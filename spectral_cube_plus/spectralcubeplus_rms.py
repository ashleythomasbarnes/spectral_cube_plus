"""Spatial RMS maps from selected channels or robust spectral clipping."""
import warnings
import numpy as np
import astropy.units as u
from astropy import stats
from astropy.io import fits
from scipy.special import erf
from scipy.optimize import minimize


def _rms_product(cube, values):
    """Package numeric noise values as a Quantity and a spatial FITS image.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Source of data unit, celestial WCS, and available common beam.
    values : numpy.ndarray
        Numeric noise map, shape (y, x), expressed in the cube data unit.

    Returns
    -------
    rms : astropy.units.Quantity
        Noise values in the cube data unit.
    rms_hdu : astropy.io.fits.PrimaryHDU
        Numeric values with celestial WCS, BUNIT, and available beam keywords.
    """
    rms = values * cube.unit
    header = cube.wcs.celestial.to_header()
    header['BUNIT'] = cube.unit.to_string(format='fits')
    if cube._meta.get('beam') is not None:
        header.update(cube.beam.to_header_keywords())
    return rms, fits.PrimaryHDU(rms.value, header)


def _mad(values, axis=None):
    # All-invalid sightlines intentionally return NaN, without a noisy warning.
    """Compute Gaussian-scaled MAD while ignoring NaNs and empty-slice warnings.

    Parameters
    ----------
    values : array-like
        Samples; invalid values should already be represented by NaN.
    axis : int or None, optional
        Reduction axis. None reduces the whole array.

    Returns
    -------
    noise : numpy.ndarray or numpy scalar
        ``astropy.stats.mad_std(..., ignore_nan=True)`` result. All-NaN
        slices produce NaN without empty-slice RuntimeWarnings.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='All-NaN slice encountered', category=RuntimeWarning)
        warnings.filterwarnings('ignore', message='Mean of empty slice', category=RuntimeWarning)
        return stats.mad_std(values, axis=axis, ignore_nan=True)


def _valid_data(cube):
    """Load cube values and replace invalid or excluded samples with NaN.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Input in (spectral, y, x) order, with an optional inclusion mask.

    Returns
    -------
    data : numpy.ndarray
        Unitless copy in the cube data unit, with nonfinite and currently
        masked samples replaced by NaN. The input cube is unchanged.
    """
    data = cube.unmasked_data[:].value
    valid = np.isfinite(data)
    if cube.mask is not None:
        valid &= cube.mask.include()
    return np.where(valid, data, np.nan)


def get_rms(cube, rms_velocities):
    """Estimate a spatial noise map from selected channel centres using MAD.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Input in (spectral, y, x) order. Existing masked and nonfinite
        samples are ignored.
    rms_velocities : array-like or astropy.units.Quantity
        One pair, shape (2,), or a collection, shape (N, 2), of finite,
        ordered bounds. Numeric bounds use the current spectral-axis unit;
        Quantity bounds must be convertible to that unit. Selection is
        inclusive at channel centres and accepts ascending or descending
        axes. Overlaps are merged so each channel is counted only once.

    Returns
    -------
    rms : astropy.units.Quantity
        Gaussian-scaled median absolute deviation along selected channels,
        shape (y, x), in the cube data unit. All-invalid sightlines are NaN.
    rms_hdu : astropy.io.fits.PrimaryHDU
        Numeric RMS image with celestial WCS, BUNIT, and available beam.

    Raises
    ------
    ValueError
        If bounds are malformed, nonfinite, reversed, incompatible with
        spectral units, or their union contains no channel centre.

    Notes
    -----
    Choose channels representative of line-free noise. Despite the name,
    this computes Gaussian-scaled MAD, not a root-mean-square about zero.
    The complete cube is materialized before selecting channels. Neither
    input data nor masks are modified and no attributes are attached.
    """
    try:
        intervals = u.Quantity(rms_velocities)
        if intervals.unit == u.dimensionless_unscaled:
            intervals = intervals.value * cube.spectral_axis.unit
        intervals = intervals.to(cube.spectral_axis.unit)
    except (TypeError, ValueError, u.UnitsError) as exc:
        raise ValueError('RMS intervals must be pairs in spectral-axis units.') from exc
    if intervals.shape == (2,):
        intervals = intervals.reshape(1, 2)
    if (intervals.ndim != 2 or intervals.shape[1] != 2 or len(intervals) == 0
            or not np.all(np.isfinite(intervals.value))
            or np.any(intervals[:, 0] > intervals[:, 1])):
        raise ValueError('RMS intervals must contain finite, ordered pairs.')
    axis = cube.spectral_axis
    selected = np.zeros(len(axis), dtype=bool)
    for lower, upper in intervals:
        selected |= (axis >= lower) & (axis <= upper)
    if not selected.any():
        raise ValueError('RMS intervals contain no cube channels.')
    data = _valid_data(cube)[selected]
    return _rms_product(cube, _mad(data, axis=0))


def erf0(x, erftarg):
    """Evaluate the absolute residual of an error-function target.

    Parameters
    ----------
    x : float or array-like
        Dimensionless argument to ``scipy.special.erf``.
    erftarg : float
        Target value of the error function.

    Returns
    -------
    residual : float or numpy.ndarray
        ``abs(erf(x) - erftarg)``. Used as the minimization objective in
        automatic noise estimation; x is not itself a Gaussian sigma count.
    """
    return abs(erf(x) - erftarg)


def get_rms_auto(cube):
    """Estimate noise from reflected negative samples and upper clipping.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Input in (spectral, y, x) order. Only finite samples included by
        its current mask contribute to noise estimates.

    Returns
    -------
    rms : astropy.units.Quantity
        Spatial noise map, shape (y, x), in the cube data unit. All-invalid
        sightlines return NaN; constant spectra may return zero.
    rms_hdu : astropy.io.fits.PrimaryHDU
        Numeric RMS image with celestial WCS, BUNIT, and available beam.

    Notes
    -----
    Start with Gaussian-scaled MAD of the finite spectrum. For each sightline
    with at least 10 negative samples, reflect those samples about zero and
    compute their MAD scale. Estimate a clipping factor by minimizing
    ``abs(erf(x) - (1 - 4 / n_channels))``, starting at x=3. Recompute MAD
    using finite samples strictly below that factor times the reflected
    noise scale. The channel count is the full cube length, not a per-pixel
    valid-sample count. Sightlines with fewer than 10 negative samples keep
    the initial estimate; an empty clipped selection produces NaN.

    This prescription assumes noise is approximately symmetric around zero
    and negative samples are noise dominated. Absorption, baseline offsets,
    and emission filling most channels can bias the result. It loads the
    full cube and loops over sightlines. Inputs are unchanged and no derived
    attributes are attached. Use explicit line-free intervals when suitable.
    """
    data = _valid_data(cube)
    channels, nx, ny = data.shape
    sig_false = minimize(erf0, 3, args=(1 - 4 / channels,)).x[0]
    rms_map = _mad(data, axis=0)
    for i in range(nx):
        for j in range(ny):
            spec = data[:, i, j]
            finite = spec[np.isfinite(spec)]
            negative = finite[finite < 0]
            if len(negative) < 10:
                continue
            sigma = _mad(np.concatenate((negative, -negative)))
            selected = finite[finite < sig_false * sigma]
            rms_map[i, j] = _mad(selected) if selected.size else np.nan
    return _rms_product(cube, rms_map)
