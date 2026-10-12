"""Independent numerical checks of requested centroid/dispersion errors."""
import numpy as np
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

from spectral_cube_plus import SpectralCubePlus
from spectral_cube_plus.spectralcubeplus_moments import get_momentmaps


def make_cube(data, step=1000, origin=0):
    data = np.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data[:, None, None]
    wcs = WCS(naxis=3)
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'VRAD']
    wcs.wcs.cunit = ['deg', 'deg', 'm/s']
    wcs.wcs.cdelt = [-0.001, 0.001, step]
    wcs.wcs.crpix = [1, 1, 1]
    wcs.wcs.crval = [0, 0, origin]
    return SpectralCubePlus(data=data*u.K, wcs=wcs)


def moments(intensity, coordinate):
    total = intensity.sum(axis=-1)
    centroid = (intensity * coordinate).sum(axis=-1) / total
    dispersion = np.sqrt((intensity * (coordinate - centroid[..., None])**2).sum(axis=-1) / total)
    return np.stack([centroid, dispersion], axis=-1)


def numerical_errors(intensity, coordinate, rms):
    # Subtract a constant origin to avoid cancellation for large frequency axes.
    coordinate = coordinate - coordinate[0]
    derivatives = []
    for i in range(len(intensity)):
        delta = np.zeros_like(intensity)
        delta[i] = 1.e-4
        derivatives.append((moments(intensity + delta, coordinate) -
                            moments(intensity - delta, coordinate)) / (2.e-4))
    return rms * np.sqrt(np.sum(np.asarray(derivatives)**2, axis=0))


@pytest.mark.parametrize('step', [1000, -1000])
@pytest.mark.parametrize('unit', [u.m/u.s, u.km/u.s, u.Hz])
def test_derivatives_units_slabs_and_headers(step, unit):
    intensity = np.array([4., 8., 17., 11., 5.])
    c = make_cube(intensity, step, origin=-2000)
    if unit == u.Hz:
        c = c.with_spectral_unit(unit, velocity_convention='radio', rest_value=100*u.GHz)
        bounds = None
        selected = c
    else:
        c = c.with_spectral_unit(unit)
        bounds = [-3, 1]*u.km/u.s
        selected = c.spectral_slab(*bounds)
    header = c.wcs.celestial.to_header()
    header['BUNIT'] = 'mK'
    rms = fits.PrimaryHDU(np.full(c.shape[1:], 500.), header)
    old_header = rms.header.copy()
    before = c.unmasked_data[:].value.copy()
    result = get_momentmaps(c, rms, bounds, velomoms=True)
    values = selected.filled_data[:, 0, 0].value
    coordinates = selected.spectral_axis.to_value(unit)
    expected = numerical_errors(values, coordinates, .5)
    for key, value in zip(['mom1err', 'sigmaerr'], expected):
        np.testing.assert_allclose(result[key].data, value, rtol=2.e-4)
        assert u.Unit(result[key].header['BUNIT']) == unit
        assert WCS(result[key].header).celestial.to_header() == c.wcs.celestial.to_header()
    for name, numerator, denominator in [('mom1s2n', 'mom1', 'mom1err'),
                                         ('sigmas2n', 'sigma', 'sigmaerr')]:
        np.testing.assert_allclose(result[name].data, result[numerator].data/result[denominator].data)
        assert u.Unit(result[name].header['BUNIT']) == u.dimensionless_unscaled
    assert rms.header == old_header
    np.testing.assert_array_equal(c.unmasked_data[:].value, before)
    assert c.mom0 is None
    disabled = get_momentmaps(c, rms)
    for key in ['mom1', 'sigma', 'mom1err', 'sigmaerr', 'mom1s2n', 'sigmas2n']:
        assert disabled[key] is None


def test_coordinate_shift_and_noise_simulation():
    intensity = np.array([40., 80., 170., 110., 50.])
    c = make_cube(intensity).with_spectral_unit(u.km/u.s)
    shifted = make_cube(intensity, origin=-10000).with_spectral_unit(u.km/u.s)
    rms = np.ones((1, 1))*u.K
    a = get_momentmaps(c, rms, velomoms=True)
    b = get_momentmaps(shifted, rms, velomoms=True)
    for key in ['mom1err', 'sigmaerr', 'sigmas2n']:
        np.testing.assert_allclose(a[key].data, b[key].data)
    np.testing.assert_allclose(b['mom1s2n'].data-a['mom1s2n'].data, -10/a['mom1err'].data)
    assert b['mom1s2n'].data.item() < 0
    rng = np.random.default_rng(342)
    samples = intensity + rng.normal(size=(50000, len(intensity)))
    measured = moments(samples, c.spectral_axis.value).std(axis=0, ddof=1)
    np.testing.assert_allclose(measured, [a['mom1err'].data.item(), a['sigmaerr'].data.item()], rtol=.02)


def test_masks_nan_signed_intensities_and_degenerate_spectra():
    # Columns: ordinary masked, NaN, empty, cancelling, one channel, all-negative,
    # mixed-sign with positive variance, and negative variance.
    data = np.array([[4, np.nan, 4, 1, 0, -4, 4, -1],
                     [8, 8, 8, -1, 5, -8, -1, 3],
                     [17, 17, 17, 0, 0, -17, 17, -1]], dtype=float)[:, None, :]
    c = make_cube(data).with_spectral_unit(u.km/u.s)
    mask = np.ones(c.shape, dtype=bool)
    mask[0, 0, 0] = False
    mask[:, 0, 2] = False
    mask[[0, 2], 0, 4] = False
    c = c.with_mask(mask)
    result = get_momentmaps(c, np.ones(c.shape[1:])*u.K, velomoms=True)
    for col in [0, 1, 5, 6]:
        valid = mask[:, 0, col] & np.isfinite(data[:, 0, col])
        expected = numerical_errors(data[valid, 0, col], c.spectral_axis.value[valid], 1)
        np.testing.assert_allclose([result[k].data[0, col] for k in ['mom1err', 'sigmaerr']], expected)
    for key in ['mom1err', 'sigmaerr', 'mom1s2n', 'sigmas2n']:
        assert np.isnan(result[key].data[0, [2, 3]]).all()
    assert result['mom1err'].data[0, 4] == 0
    for key in ['sigmaerr', 'mom1s2n', 'sigmas2n']:
        assert np.isnan(result[key].data[0, 4])
    assert np.isfinite(result['mom1err'].data[0, 7])
    assert np.isnan(result['sigmaerr'].data[0, 7])
    single = make_cube([5.], origin=2000)
    single_result = get_momentmaps(single, np.ones((1, 1))*u.K, velomoms=True)
    assert single_result['mom1err'].data.item() == 0
    assert np.isnan(single_result['sigmaerr'].data.item())
    # Zero-valued included channels still contribute noise derivatives.
    zero_width = get_momentmaps(make_cube([0, 5, 0]), np.ones((1, 1))*u.K,
                               velomoms=True)
    assert zero_width['mom1err'].data.item() > 0
    assert np.isnan(zero_width['sigmaerr'].data.item())
    # A non-NaN fill value must not add excluded channels to the error sums.
    filled = get_momentmaps(c.with_fill_value(0), np.ones(c.shape[1:])*u.K,
                           velomoms=True)
    for key in ['mom1err', 'sigmaerr']:
        np.testing.assert_allclose(filled[key].data[0, 0], result[key].data[0, 0])


def test_invalid_rms():
    c = make_cube([4, 8, 17])
    for rms in [0, np.nan, np.inf]:
        result = get_momentmaps(c, np.full((1, 1), rms)*u.K, velomoms=True)
        for key in ['mom1err', 'sigmaerr', 'mom1s2n', 'sigmas2n']:
            assert np.isnan(result[key].data).all()
    with pytest.raises(ValueError, match='nonnegative'):
        get_momentmaps(c, -np.ones((1, 1))*u.K, velomoms=True)
    with pytest.raises(ValueError, match='shape'):
        get_momentmaps(c, np.ones((2, 1))*u.K, velomoms=True)
    with pytest.raises(u.UnitsError):
        get_momentmaps(c, np.ones((1, 1))*u.Jy, velomoms=True)
