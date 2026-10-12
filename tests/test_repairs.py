import numpy as np
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import Projection

from spectral_cube_plus import SpectralCubePlus
from spectral_cube_plus.spectralcubeplus_moments import get_momentmaps
from spectral_cube_plus.spectralcubeplus_masks import (
    get_prunemask, get_prunemaskvelo, get_expmask, _beam_geometry,
)


def cube(data=None, step=1000):
    wcs = WCS(naxis=3)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'VRAD']
    wcs.wcs.cunit = ['deg', 'deg', 'm/s']
    wcs.wcs.cdelt = [-0.001, 0.002, step]
    wcs.wcs.crpix = [1, 1, 1]
    if data is None:
        data = np.ones((12, 3, 3))
    header = wcs.to_header()
    header['BMAJ'] = 0.003
    header['BMIN'] = 0.002
    header['BPA'] = 0
    header['OBSERVER'] = 'regression'
    return SpectralCubePlus(data=data * u.K, wcs=wcs, header=header)


def rms_hdu(c, value=2):
    header = c.wcs.celestial.to_header()
    header['BUNIT'] = 'K'
    return fits.PrimaryHDU(np.full(c.shape[1:], value, dtype=float), header)


@pytest.mark.parametrize('kind', ['quantity', 'hdu', 'projection'])
def test_error_units_masks_and_repeat(kind, tmp_path):
    data = np.ones((12, 3, 3))
    data[0, 0, 0] = np.nan
    c = cube(data).with_spectral_unit(u.km / u.s)
    mask = np.ones(c.shape, bool)
    mask[:4] = False
    mask[:, 2, 2] = False
    c = c.with_mask(mask)
    hdu = rms_hdu(c)
    original_header = hdu.header.copy()
    rms = {'quantity': 2000 * u.mK * np.ones((3, 3)),
           'hdu': hdu, 'projection': Projection.from_hdu(hdu)}[kind]
    for _ in range(2):
        error = c.moment0err(rms=rms)
        assert isinstance(error, Projection)
        assert error.unit == u.K * u.km / u.s
        np.testing.assert_allclose(error.value[:2], 2 * np.sqrt(8))
        assert np.isnan(error.value[2, 2])
        error.write(tmp_path / 'error.fits', overwrite=True)
    assert hdu.header == original_header
    assert callable(c.moment0err)


def test_momentmaps_slab_and_snr():
    c = cube()
    result = get_momentmaps(c, rms_hdu(c), [0, 2], velomoms=True)
    np.testing.assert_allclose(result['mom0'].data, 3000)
    np.testing.assert_allclose(result['mom0err'].data, 2000 * np.sqrt(3))
    assert u.Unit(result['mom0s2n'].header.get('BUNIT', '')) == u.dimensionless_unscaled
    assert result['mom1'] is not None
    converted = c.with_spectral_unit(u.km/u.s)
    other = get_momentmaps(converted, rms_hdu(converted), [0*u.m/u.s, 2000*u.m/u.s])
    np.testing.assert_allclose(other['mom0err'].data * 1000, result['mom0err'].data)
    with pytest.raises(ValueError):
        get_momentmaps(c, rms_hdu(c), [2, 0])


def test_rms_nan_fallback_and_write(tmp_path):
    c = cube(np.array([-2., -1., 0., 1., 2., np.nan])[:, None, None])
    np.testing.assert_allclose(c.get_rms_auto().value, 1.4826022185)
    c.rms_hdu.writeto(tmp_path / 'rms.fits')
    with fits.open(tmp_path / 'rms.fits') as hdus:
        hdu = hdus[0]
        assert hdu.data.shape == (1, 1)
        assert u.Unit(hdu.header['BUNIT']) == u.K
        assert WCS(hdu.header).pixel_n_dim == 2
    empty = cube(np.full((12, 1, 1), np.nan))
    assert np.isnan(empty.get_rms_auto().value).all()


def test_rms_intervals_union():
    c = cube(np.arange(12.)[:, None, None])
    overlap = c.get_rms(rms_velocities=[[0, 4000], [2000, 6000]])
    full = c.get_rms(rms_velocities=[0, 6000])
    np.testing.assert_allclose(overlap, full)
    quantity = c.get_rms(rms_velocities=np.array([[0, 6]]) * u.km/u.s)
    np.testing.assert_allclose(quantity, full)
    with pytest.raises(ValueError):
        c.get_rms(rms_velocities=[5, 1])


def test_pruning_labels_and_edges():
    mask = np.zeros((1, 9, 9), bool)
    mask[0, 1, 1:8] = True
    mask[0, 1:8, 1] = True
    mask[0, 3:7, 3:7] = True
    assert get_prunemask(mask, 15).sum() == 16
    edge = np.array([1, 1, 0, 0, 1], bool)[:, None, None]
    assert not get_prunemaskvelo(edge).any()
    assert get_prunemaskvelo(np.ones((5, 1, 1), bool), 2).sum() == 1
    with pytest.raises(ValueError):
        get_prunemaskvelo(edge, -1)


def test_geometry_and_existing_mask():
    # Gaussian area is 3.4 pixels; the old formula gives 9.4 pixels.
    c = cube(np.full((15, 15, 15), 10.))
    mask = np.ones(c.shape, bool)
    mask[:, 7, 7] = False
    c = c.with_mask(mask)
    result = get_expmask(c, np.ones((15, 15))*u.K, beamarea=1)
    assert result.any()
    assert not result[:, 7, 7].any()
    area, structure = _beam_geometry(c, 1)
    np.testing.assert_allclose(area, np.pi*0.003*0.002/(4*np.log(2)*0.001*0.002))
    assert structure.shape == (5, 3, 5)
    _, large_structure = _beam_geometry(c, 10)
    assert large_structure.shape[0] > 11


def test_gaussian_beam_pruning_threshold():
    data = np.zeros((5, 9, 9))
    data[:, 3:5, 3:5] = 10
    c = cube(data)
    # Four pixels survive one Gaussian beam (3.4 pixels), but not two.
    one = get_expmask(c, np.ones((9, 9))*u.K, beamarea=1, radfrac=0, npix=0)
    two = get_expmask(c, np.ones((9, 9))*u.K, beamarea=2, radfrac=0, npix=0)
    assert one.sum() == 20
    assert two.sum() == 0


def test_moment_options_and_reset():
    c = cube()
    c.moment0(how='cube')
    cached = c.mom0
    c.moment0(axis=1, how='slice')
    assert c.mom0 is cached
    c.get_rms_auto()
    changed = c.with_spectral_unit(u.km/u.s)
    assert changed.rms is None and changed.mom0 is None
    assert c.with_mask(np.ones(c.shape, bool)).rms is None


@pytest.mark.parametrize('temperature_unit', ['K (Tmb)', 'K (Ta*)'])
def test_read_metadata_units_and_beam(tmp_path, temperature_unit):
    c = cube()
    header = c.header
    header['BUNIT'] = temperature_unit
    header['RESTFRQ'] = 1e11
    path = tmp_path / 'input.fits'
    fits.PrimaryHDU(c._data, header).writeto(path)
    for _ in range(2):
        read = SpectralCubePlus.read(path)
        assert read.unit.is_equivalent(u.K)
        assert read.header['OBSERVER'] == 'regression'
        assert read.header['RESTFRQ'] == 1e11
        assert read.beam == c.beam
        np.testing.assert_allclose(read.wcs.wcs.cdelt, c.wcs.wcs.cdelt)
        np.testing.assert_allclose(read.spectral_axis.to_value(u.m/u.s), c.spectral_axis.value)
        read.write(tmp_path / 'roundtrip.fits', overwrite=True)
        assert fits.getheader(tmp_path / 'roundtrip.fits')['OBSERVER'] == 'regression'


def test_bad_rms_alignment_and_missing_rms():
    c = cube()
    with pytest.raises(ValueError):
        c.moment0err()
    with pytest.raises(ValueError):
        c.moment0err(rms=np.ones((2, 2))*u.K)
    with pytest.raises(u.UnitsError):
        c.moment0err(rms=np.ones((3, 3))*u.Jy)
    hdu = rms_hdu(c)
    hdu.header['CRVAL1'] = 1
    with pytest.raises(ValueError):
        c.moment0err(rms=hdu)
    with pytest.raises(ValueError):
        c.moment0err(rms=-np.ones((3, 3))*u.K)


def test_invalid_rms_and_mask_prerequisites():
    c = cube()
    rms = np.ones((3, 3))*u.K
    rms[0, 0] = np.inf*u.K
    assert np.isnan(c.moment0err(rms=rms).value[0, 0])
    assert not c.get_expmask().any()
    assert c.rms is not None
    no_beam = SpectralCubePlus(np.ones((12, 3, 3))*u.K, wcs=c.wcs)
    with pytest.raises(ValueError, match='beam'):
        get_expmask(no_beam, rms)
    with pytest.raises(ValueError):
        get_expmask(c, rms, hthresh=1, lthresh=2)


def test_finite_fill_value_and_invalid_data():
    data = np.ones((12, 3, 3))
    data[0, 0, 0] = np.nan
    data[1, 0, 0] = np.inf
    c = cube(data).with_fill_value(0)
    mask = np.ones(c.shape, bool)
    mask[-2:] = False
    c = c.with_mask(mask)
    error = c.moment0err(rms=np.ones((3, 3))*u.K)
    np.testing.assert_allclose(error.value[0, 0], 1000*np.sqrt(8))
    np.testing.assert_allclose(error.value[1, 1], 1000*np.sqrt(10))
    np.testing.assert_allclose(c.get_rms_auto().value, 0)


def test_nonuniform_spacing_rejected():
    # A radio velocity conversion of a linear wavelength axis is nonlinear.
    wcs = WCS(naxis=3)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'WAVE']
    wcs.wcs.cunit = ['deg', 'deg', 'm']
    wcs.wcs.crpix = [1, 1, 1]
    wcs.wcs.crval = [0, 0, 0.001]
    wcs.wcs.cdelt = [-0.001, 0.001, 0.0001]
    c = SpectralCubePlus(np.ones((12, 3, 3))*u.K, wcs=wcs)
    c = c.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=0.001*u.m)
    with pytest.raises(ValueError, match='uniform'):
        c.moment0err(rms=np.ones((3, 3))*u.K)


def test_negative_channel_step_and_single_channel():
    c = cube(step=-1000)
    error = c.moment0err(rms=np.ones((3, 3))*u.K)
    np.testing.assert_allclose(error.value, 1000*np.sqrt(12))
    single = c.spectral_slab(0*u.m/u.s, 0*u.m/u.s)
    np.testing.assert_allclose(single.moment0err(rms=np.ones((3, 3))*u.K).value, 1000)


def test_expanded_cube_receives_rms():
    c = cube(np.random.default_rng(3).normal(size=(30, 15, 15)))
    c.get_rms_auto()
    c.get_expmask()
    assert c.masked.rms is c.rms
    before = c.rms_hdu.header.copy()
    c.masked.moment0err()
    assert c.rms_hdu.header == before


def test_read_rejects_per_channel_beams(tmp_path):
    from spectral_cube import VaryingResolutionSpectralCube
    from radio_beam import Beams
    c = cube()
    varying = VaryingResolutionSpectralCube(
        c._data*u.K, wcs=c.wcs,
        beams=Beams(major=np.full(12, 3)*u.arcsec,
                    minor=np.full(12, 2)*u.arcsec, pa=np.zeros(12)*u.deg))
    path = tmp_path / 'beams.fits'
    varying.write(path)
    with pytest.raises(NotImplementedError, match='Per-channel'):
        SpectralCubePlus.read(path)


def test_read_dask_backend(tmp_path):
    from spectral_cube import DaskSpectralCube
    path = tmp_path / 'dask.fits'
    cube().write(path)
    c = SpectralCubePlus.read(path, use_dask=True)
    assert isinstance(c, DaskSpectralCube)
    assert isinstance(c, SpectralCubePlus)
    np.testing.assert_allclose(c.moment0().value, 12000)
    c.get_rms_auto()
    assert c.with_spectral_unit(u.km/u.s).rms is None
