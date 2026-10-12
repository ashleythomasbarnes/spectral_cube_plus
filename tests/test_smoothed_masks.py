import numpy as np
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from radio_beam import Beam
from spectral_cube import DaskSpectralCube, Projection

from spectral_cube_plus import SpectralCubePlus
from spectral_cube_plus.spectralcubeplus_masks import get_expmask
from spectral_cube_plus.spectralcubeplus_rms import get_rms_auto


def emission_cube(unit=u.K):
    wcs = WCS(naxis=3)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'VRAD']
    wcs.wcs.cunit = ['deg', 'deg', 'm/s']
    wcs.wcs.cdelt = [-1/3600, 1/3600, 1000]
    wcs.wcs.crpix = [18, 18, 1]
    data = np.random.default_rng(12).normal(size=(24, 35, 35))
    y, x = np.mgrid[:35, :35]
    data[8:14] += 8 * np.exp(-((x-17)**2 + (y-17)**2)/50)
    data[10, 17, 17] += 40
    return SpectralCubePlus(data * unit, wcs=wcs,
                            beam=Beam(major=3*u.arcsec, minor=2*u.arcsec,
                                      pa=25*u.deg)).with_spectral_unit(u.km/u.s)


OPTIONS = dict(hthresh=3, lthresh=2, beamarea=0, radfrac=0, npix=1)


def test_auto_and_stored_rms(capsys, monkeypatch):
    c = emission_cube()
    assert c.mask_cube is None
    mask = c.get_expmask(**OPTIONS)
    expected, _ = get_rms_auto(c)
    np.testing.assert_allclose(c.rms.value, expected.value)
    np.testing.assert_array_equal(mask, get_expmask(c, expected, **OPTIONS))
    output = capsys.readouterr().out
    assert 'calculating native RMS automatically' in output
    assert 'median=' in output and 'range=[' in output and 'K' in output
    stored = c.rms
    def unexpected_auto(*args, **kwargs):
        pytest.fail('Stored RMS should be reused')
    monkeypatch.setattr('spectral_cube_plus.spectral_cube_plus._get_rms_auto', unexpected_auto)
    # Existing positional threshold parameters remain supported.
    c.get_expmask(3, 2, 0, 0, 1)
    assert c.rms is stored and c.masked.rms is stored
    assert 'Using stored native RMS' in capsys.readouterr().out
    assert c.mask_cube.beam == c.beam
    assert c.with_mask(mask).mask_cube is None
    assert c.spectral_slab(8*u.km/u.s, 13*u.km/u.s).mask_cube is None


@pytest.mark.parametrize('kind', ['quantity', 'hdu', 'projection'])
def test_explicit_rms(kind, capsys):
    c = emission_cube()
    c.get_rms_auto()
    header = c.wcs.celestial.to_header()
    header['BUNIT'] = 'mK'
    hdu = fits.PrimaryHDU(np.full(c.shape[1:], 2000.), header)
    original = hdu.header.copy()
    rms = {'quantity': hdu.data*u.mK, 'hdu': hdu,
           'projection': Projection.from_hdu(hdu)}[kind]
    c.get_expmask(rms=rms, **OPTIONS)
    np.testing.assert_allclose(c.rms.value, 2)
    assert c.rms.unit == u.K
    assert c.masked.rms is c.rms
    assert c.masked.rms_hdu is not c.rms_hdu
    np.testing.assert_allclose(c.mask_cube.rms.value, 2)
    np.testing.assert_allclose(c.rms_hdu.data, 2)
    assert u.Unit(c.rms_hdu.header['BUNIT']) == u.K
    assert hdu.header == original
    assert 'Using explicit native RMS: median=2' in capsys.readouterr().out


@pytest.mark.parametrize('unit', [u.K, u.Jy/u.beam])
def test_smoothed_detection_and_native_moments(unit, capsys):
    c = emission_cube(unit)
    c._data[9, 17, 16] = np.nan
    c._data[10, 17, 18] = np.inf
    included = np.ones(c.shape, bool)
    included[:, 16, 17] = False
    c = c.with_mask(included)
    # An interval estimate must remain the native uncertainty reference.
    native = c.get_rms([0, 6]*u.km/u.s)
    native_hdu = c.rms_hdu
    original_data = c._data.copy()
    original_mask = c.mask.include().copy()
    mask = c.get_expmask(smooth_factor=2, **OPTIONS)
    detection = c.mask_cube
    assert detection.beam == Beam(major=6*u.arcsec, minor=4*u.arcsec, pa=c.beam.pa)
    assert detection.shape == c.shape
    assert detection.wcs.wcs.compare(c.wcs.wcs)
    np.testing.assert_allclose(detection.spectral_axis, c.spectral_axis)
    assert c.rms is native and c.rms_hdu is native_hdu
    assert c.masked.rms is native
    assert c.masked.beam == c.beam and c.masked.unit == unit
    np.testing.assert_array_equal(c._data, original_data)
    np.testing.assert_array_equal(c.mask.include(), original_mask)
    np.testing.assert_array_equal(c.masked._data, original_data)
    assert not mask[:, 16, 17].any()
    assert not mask[9, 17, 16] and not mask[10, 17, 18]
    assert mask.any()
    valid = included & np.isfinite(original_data)
    manual = c.with_mask(valid).convolve_to(detection.beam)
    manual_rms, _ = get_rms_auto(manual)
    np.testing.assert_allclose(detection.rms.value, manual_rms.value, equal_nan=True)
    np.testing.assert_array_equal(mask, get_expmask(manual, manual_rms, **OPTIONS) & valid)
    assert detection.rms is not native
    if unit == u.K:
        assert np.nanmedian(detection.rms.value) < np.nanmedian(native.value)
    assert detection.rms_hdu.header['BMAJ'] == pytest.approx(6/3600)
    assert c.masked.rms_hdu.header['BMAJ'] == pytest.approx(3/3600)
    expected_mom0 = np.where(mask, original_data, 0).sum(axis=0)
    counts = mask.sum(axis=0)
    expected_mom0[counts == 0] = np.nan
    np.testing.assert_allclose(c.masked.moment0().value, expected_mom0, equal_nan=True)
    expected_error = native.value * np.sqrt(counts)
    expected_error[counts == 0] = np.nan
    np.testing.assert_allclose(c.masked.moment0err().value, expected_error, equal_nan=True)
    output = capsys.readouterr().out
    assert 'Using stored native RMS' in output
    assert 'Retaining native RMS for uncertainties' in output
    assert 'Smoothed detection RMS: median=' in output
    assert c.masked.with_spectral_unit(u.m/u.s).mask_cube is None


def test_smoothing_without_stored_rms(capsys):
    c = emission_cube()
    assert c.mask is None
    mask = c.get_expmask(smooth_factor=2, **OPTIONS)
    assert mask.any()
    native, _ = get_rms_auto(c)
    np.testing.assert_allclose(c.rms.value, native.value)
    assert c.masked.rms is c.rms
    assert 'calculating native RMS automatically' in capsys.readouterr().out


def test_factor_one_recomputes_detection_only(monkeypatch):
    c = emission_cube()
    c.get_rms([0, 6]*u.km/u.s)
    stored = c.rms
    def unexpected_convolution(*args, **kwargs):
        pytest.fail('Factor one must skip convolution')
    monkeypatch.setattr(SpectralCubePlus, 'convolve_to', unexpected_convolution)
    c.get_expmask(smooth_factor=1, **OPTIONS)
    expected, _ = get_rms_auto(c)
    assert c.rms is stored and c.mask_cube.rms is not stored
    assert c.mask_cube.beam == c.beam
    np.testing.assert_allclose(c.mask_cube.rms.value, expected.value)


def test_masked_fill_values_do_not_leak_into_detection():
    c = emission_cube()
    included = np.ones(c.shape, bool)
    included[:, 15:19, 15:19] = False
    c = c.with_mask(included)
    filled = c.with_fill_value(1e6)
    options = dict(hthresh=3, lthresh=2, beamarea=1, radfrac=0.3, npix=1)
    mask = filled.get_expmask(smooth_factor=2, **options)
    expected = c.get_expmask(smooth_factor=2, **options)
    np.testing.assert_array_equal(mask, expected)
    np.testing.assert_allclose(filled.mask_cube.rms.value, c.mask_cube.rms.value,
                               equal_nan=True)
    assert filled.masked._fill_value == 1e6
    # Beam-scaled closing and pruning must use the detection beam.
    np.testing.assert_array_equal(mask, get_expmask(
        filled.mask_cube, filled.mask_cube.rms, **options) & included)


@pytest.mark.parametrize('kind', ['shape', 'unit', 'wcs'])
def test_bad_explicit_rms_does_not_replace_products(kind):
    c = emission_cube()
    c.get_expmask(**OPTIONS)
    previous = (c.rms, c.rms_hdu, c.mask_cube, c.masked, c.expmask)
    if kind == 'shape':
        rms = np.ones((2, 2))*u.K
    elif kind == 'unit':
        rms = np.ones(c.shape[1:])*u.Jy
    else:
        rms = c.rms_hdu.copy()
        rms.header['CRVAL1'] = 1
    with pytest.raises((ValueError, u.UnitsError)):
        c.get_expmask(rms=rms, **OPTIONS)
    for current, old in zip((c.rms, c.rms_hdu, c.mask_cube, c.masked, c.expmask), previous):
        assert current is old


@pytest.mark.parametrize('factor', [0, -1, 0.5, np.nan, np.inf, True, np.bool_(False), '2', 2j])
def test_invalid_factors_leave_state_unchanged(factor, capsys):
    c = emission_cube()
    with pytest.raises(ValueError, match='smooth_factor'):
        c.get_expmask(smooth_factor=factor)
    assert c.rms is None and c.mask_cube is None and c.masked is None
    assert capsys.readouterr().out == ''


def test_conflict_and_missing_beam_leave_state_unchanged(capsys):
    c = emission_cube()
    c.get_expmask(**OPTIONS)
    previous = (c.rms, c.rms_hdu, c.mask_cube, c.masked, c.expmask)
    capsys.readouterr()
    with pytest.raises(ValueError, match='cannot both'):
        c.get_expmask(rms=c.rms, smooth_factor=2)
    for current, old in zip((c.rms, c.rms_hdu, c.mask_cube, c.masked, c.expmask), previous):
        assert current is old
    assert capsys.readouterr().out == ''
    no_beam = SpectralCubePlus(c._data*u.K, wcs=c.wcs)
    with pytest.raises(ValueError, match='beam'):
        no_beam.get_expmask(smooth_factor=2)
    assert no_beam.rms is None and no_beam.mask_cube is None


def test_all_invalid_rms_summary(capsys):
    c = emission_cube()
    c = c.with_mask(np.zeros(c.shape, bool))
    assert not c.get_expmask(**OPTIONS).any()
    assert 'no finite values (K)' in capsys.readouterr().out


def test_dask_smoothed_mask(tmp_path):
    c = emission_cube()
    path = tmp_path/'cube.fits'
    c.write(path)
    dask = SpectralCubePlus.read(path, use_dask=True)
    mask = dask.get_expmask(smooth_factor=2, **OPTIONS)
    assert mask.any()
    assert isinstance(dask.mask_cube, DaskSpectralCube)
    assert isinstance(dask.masked, DaskSpectralCube)
    assert dask.masked.beam == dask.beam
    np.testing.assert_array_equal(mask, c.get_expmask(smooth_factor=2, **OPTIONS))
    np.testing.assert_allclose(dask.masked.moment0().value, c.masked.moment0().value,
                               equal_nan=True)
    np.testing.assert_allclose(dask.masked.moment0err().value, c.masked.moment0err().value,
                               equal_nan=True)
