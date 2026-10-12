"""Moment products and independent-channel uncertainties."""
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import Projection, wcs_utils


def _rms_quantity(rms, unit, shape, wcs=None):
    """Validate, align, and convert a spatial RMS map without modifying it.

    Parameters
    ----------
    rms : astropy.units.Quantity, spectral_cube.Projection, or FITS image HDU
        Map to normalize. FITS PrimaryHDU/ImageHDU inputs must declare BUNIT.
    unit : astropy.units.Unit
        Required output unit, convertible from the map unit.
    shape : tuple of int
        Required exact spatial shape; no broadcasting is performed.
    wcs : astropy.wcs.WCS, optional
        Target WCS. If the input has celestial WCS, compare it with the
        target celestial WCS using tolerance 1e-8.

    Returns
    -------
    values : astropy.units.Quantity
        Map in the requested unit, with nonfinite values replaced by NaN.
        Zero is allowed; finite negative values are rejected.

    Raises
    ------
    TypeError
        If the map representation is unsupported.
    ValueError
        If BUNIT is missing, the shape/WCS is misaligned, or finite RMS is
        negative.
    astropy.units.UnitsError
        If units are incompatible.
    """
    rms_wcs = None
    if isinstance(rms, Projection):
        values = u.Quantity(rms)
        rms_wcs = rms.wcs
    elif isinstance(rms, (fits.PrimaryHDU, fits.ImageHDU)):
        if 'BUNIT' not in rms.header:
            raise ValueError('RMS FITS HDU must declare BUNIT.')
        values = u.Quantity(rms.data, u.Unit(rms.header['BUNIT']))
        candidate = WCS(rms.header)
        if candidate.has_celestial:
            rms_wcs = candidate.celestial
    elif isinstance(rms, u.Quantity):
        values = rms
    else:
        raise TypeError('RMS must be a Quantity, Projection, or FITS image HDU.')
    if values.shape != shape:
        raise ValueError('RMS map shape must match the spatial cube shape.')
    values = values.to(unit)
    if np.any(values.value[np.isfinite(values.value)] < 0):
        raise ValueError('RMS values must be nonnegative.')
    if rms_wcs is not None and wcs is not None:
        if not wcs_utils.check_equality(rms_wcs.celestial, wcs.celestial,
                                        wcs_tolerance=1e-8):
            raise ValueError('RMS map WCS does not match the cube spatial WCS.')
    return u.Quantity(np.where(np.isfinite(values.value), values.value, np.nan), unit)


def _channel_width(cube):
    """Return a positive spectral pixel width consistent with cube moments.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Cube with uniformly spaced spectral coordinates. Both ascending and
        descending axes, including a single-channel axis, are supported.

    Returns
    -------
    width : astropy.units.Quantity
        Scalar channel width in the current spectral-axis unit, using the
        upstream moment pixel-width convention and unit scaling.

    Raises
    ------
    ValueError
        If spacings differ beyond rtol=1e-6, or width is nonfinite/nonpositive.
    """
    axis = cube.spectral_axis
    if len(axis) > 1:
        differences = np.diff(axis.to_value(axis.unit))
        if not np.allclose(differences, differences[0], rtol=1e-6, atol=0):
            raise ValueError('Moment errors require uniform spectral spacing.')
    # Match SpectralCube.moment's pixel-width convention, including unit scaling.
    width = cube._pix_size_slice(0) * axis.unit
    if not np.isfinite(width.value) or width.value <= 0:
        raise ValueError('Spectral channel width must be finite and positive.')
    return width


def get_nchancube(cube):
    """Count finite, included samples along each cube sightline.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Cube in (spectral, y, x) order with an optional inclusion mask.

    Returns
    -------
    nchan : numpy.ndarray of int
        Shape (y, x). Counts finite data values allowed by the current mask;
        empty sightlines have zero counts.

    Notes
    -----
    Materializes the full unmasked data to check finiteness, then applies
    mask inclusion. The cube and its mask are unchanged.
    """
    valid = np.isfinite(cube.unmasked_data[:].value)
    if cube.mask is not None:
        valid &= cube.mask.include()
    return np.sum(valid, axis=0)


def get_nchanmask(mask):
    """Count included channels in a SpectralCube mask.

    Parameters
    ----------
    mask : spectral_cube.MaskBase
        Mask whose ``include()`` gives a (spectral, y, x) boolean array.

    Returns
    -------
    nchan : numpy.ndarray of int
        Shape (y, x), counting True samples along axis 0.

    Notes
    -----
    Counts mask inclusion only; does not independently inspect data for
    nonfinite values. The input is unchanged.
    """
    return np.sum(mask.include(), axis=0)


def get_nchanbool(data):
    """Count True samples along axis zero of an inclusion array.

    Parameters
    ----------
    data : array-like
        Boolean inclusion array, normally shape (spectral, y, x). Values
        are cast to bool; pass a mask rather than numeric science data.

    Returns
    -------
    nchan : numpy.ndarray of int
        Counts over axis 0, normally shape (y, x). All-False sightlines
        return zero. The input array is unchanged.
    """
    return np.sum(np.asarray(data, dtype=bool), axis=0)


def get_mom0err(mom0, nchan, rms_map, velo_res):
    """Build a FITS uncertainty map for spectral integrated intensity.

    Parameters
    ----------
    mom0 : spectral_cube.Projection
        Moment-zero map supplying the output shape, WCS, header, and unit.
    nchan : array-like
        Nonnegative included-channel counts with exactly the moment shape.
        Count finite samples contributing to the same integration as mom0.
    rms_map : astropy.units.Quantity, spectral_cube.Projection, or FITS image HDU
        Spatial RMS map with the moment shape and units convertible to
        ``mom0.unit / velo_res.unit``. FITS inputs must declare BUNIT;
        available celestial WCS must match mom0.
    velo_res : astropy.units.Quantity
        Finite nonzero scalar spectral channel width. Its absolute value
        is used, so descending spectral axes are supported. The unit need
        not be velocity, but must be consistent with the moment unit.

    Returns
    -------
    error_hdu : astropy.io.fits.PrimaryHDU
        Numeric ``RMS * abs(velo_res) * sqrt(nchan)`` in the moment unit,
        using a copied moment header. Zero counts, nonfinite RMS, and zero
        RMS give NaN errors.

    Raises
    ------
    ValueError
        If channel width is invalid, counts have the wrong shape or negative
        values, or RMS metadata, values, or alignment are invalid.
    TypeError
        If RMS is an unsupported representation.
    astropy.units.UnitsError
        If RMS/channel-width units are incompatible with the moment unit.

    Notes
    -----
    Assumes independent channels, constant spectral width, and constant RMS
    along a sightline. Correlated channels and calibration uncertainties
    are not included. Inputs and their headers are not modified. This
    low-level helper does not verify that counts describe the supplied map.
    """
    width = np.abs(u.Quantity(velo_res))
    if width.ndim != 0 or not np.isfinite(width.value) or width.value <= 0:
        raise ValueError('Channel width must be a finite positive scalar.')
    rms = _rms_quantity(rms_map, mom0.unit / width.unit, mom0.shape, mom0.wcs)
    nchan = np.asarray(nchan)
    if nchan.shape != mom0.shape or np.any(nchan < 0):
        raise ValueError('Channel counts must be nonnegative and match moment shape.')
    error = (rms * width * np.sqrt(nchan)).to(mom0.unit).value
    error = np.where((nchan > 0) & (rms.value > 0), error, np.nan)
    header = mom0.hdu.header.copy()
    header['BUNIT'] = mom0.unit.to_string(format='fits')
    return fits.PrimaryHDU(error, header=header)


def _moment0_error(cube, rms, mom0=None):
    """Combine channel counting, width validation, and error propagation.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Cube whose finite, included samples define the integration.
    rms : astropy.units.Quantity, spectral_cube.Projection, or FITS image HDU
        Explicit aligned RMS map; None is rejected.
    mom0 : spectral_cube.Projection, optional
        Matching moment-zero map. If omitted, compute it from this cube.

    Returns
    -------
    error : spectral_cube.Projection
        Independent-channel spectral moment-zero uncertainty. Uses
        ``get_mom0err`` for units, RMS validation, and empty-sightline handling.

    Raises
    ------
    ValueError
        If RMS is missing or spectral spacing is nonuniform/invalid.
    """
    if rms is None:
        raise ValueError('Compute RMS first or supply an RMS map explicitly.')
    width = _channel_width(cube)
    if mom0 is None:
        mom0 = cube.moment(order=0, axis=0)
    return Projection.from_hdu(get_mom0err(mom0, get_nchancube(cube), rms, width))


def _velocity_errors(cube, rms, mom1, sigma):
    """Propagate channel noise through the centroid and dispersion derivatives.

    RMS is already validated and in cube units. Accumulate only spatial
    arrays, with finite included channels defining the fixed integration mask.
    """
    total = np.zeros(mom1.shape)
    count = np.zeros(mom1.shape, dtype=int)
    centroid_sum = np.zeros(mom1.shape)
    variance_sum = np.zeros(mom1.shape)
    dispersion2 = sigma.value ** 2
    for index, coordinate in enumerate(cube.spectral_axis.to_value(mom1.unit)):
        view = (index, slice(None), slice(None))
        values = cube.unmasked_data[view].value
        valid = np.isfinite(values)
        if cube.mask is not None:
            valid &= cube.mask.include(view=view)
        count += valid
        total += np.where(valid, values, 0)
        offset2 = (coordinate - mom1.value) ** 2
        centroid_sum += np.where(valid, offset2, 0)
        variance_sum += np.where(valid, (offset2 - dispersion2) ** 2, 0)
    usable = (np.isfinite(total) & (total != 0) &
              np.isfinite(rms.value) & (rms.value > 0))
    centroid_error = np.full(mom1.shape, np.nan)
    sigma_error = np.full(sigma.shape, np.nan)
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        np.divide(rms.value * np.sqrt(centroid_sum), np.abs(total),
                  out=centroid_error, where=usable & np.isfinite(mom1.value))
        np.divide(rms.value * np.sqrt(variance_sum),
                  2 * sigma.value * np.abs(total), out=sigma_error,
                  where=usable & np.isfinite(mom1.value) &
                  np.isfinite(sigma.value) & (sigma.value > 0) & (count > 1))
    # A single coordinate has an exactly fixed centroid, despite rounding in
    # the upstream weighted mean; its dispersion derivative is undefined.
    centroid_error[usable & np.isfinite(mom1.value) & (count == 1)] = 0
    outputs = []
    for values, moment in [(centroid_error, mom1), (sigma_error, sigma)]:
        header = moment.hdu.header.copy()
        header['BUNIT'] = moment.unit.to_string(format='fits')
        outputs.append(fits.PrimaryHDU(np.where(np.isfinite(values), values, np.nan), header))
    return outputs


def get_momentmaps(cube, rms, mom_velocity=None, velomoms=False):
    """Calculate FITS moments, noise products, and signal-to-noise maps.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Input in (spectral, y, x) order; respects its current inclusion mask.
        Supply the already masked cube for emission-mask products.
    rms : astropy.units.Quantity, spectral_cube.Projection, or FITS image HDU
        Spatial noise map, shape (y, x), convertible to cube data units.
        FITS inputs require BUNIT; available celestial WCS must match.
    mom_velocity : array-like or astropy.units.Quantity, optional
        Two finite ordered velocity bounds. Numeric bounds mean km/s;
        Quantity bounds must have velocity units. The cube must already
        have a velocity axis when bounds are supplied. Uses the upstream
        ``spectral_slab`` channel selection. None or the legacy empty string
        uses the full cube. A range entirely outside the cube is rejected.
    velomoms : bool, default False
        Also compute moment one, linewidth standard deviation, their
        first-order statistical uncertainties, and signed moment/error ratios.

    Returns
    -------
    products : dict of str to astropy.io.fits.PrimaryHDU or None
        Spatial maps of shape (y, x), with celestial WCS and BUNIT:

        * ``mom0``: integrated intensity, data unit times spectral unit.
        * ``max``: maximum included intensity, data unit.
        * ``mom0err``: independent-channel uncertainty in mom0 units.
        * ``maxerr``: supplied spatial RMS in the cube data unit; a noise
          reference, not an uncertainty model for the maximum statistic.
        * ``mom0s2n``: mom0 / mom0err, dimensionless.
        * ``maxs2n``: max / RMS, dimensionless.
        * ``mom1``: intensity-weighted coordinate, spectral unit; None
          unless velomoms is True.
        * ``sigma``: linewidth standard deviation, spectral unit; None
          unless velomoms is True.
        * ``mom1err``, ``sigmaerr``: centroid and dispersion uncertainties,
          spectral unit; None unless velomoms is True.
        * ``mom1s2n``, ``sigmas2n``: signed moment/error ratios, dimensionless;
          None unless velomoms is True.

    Raises
    ------
    ValueError
        If velocity bounds are invalid or outside the cube, its axis is not
        velocity when bounds are supplied, spectral spacing is nonuniform,
        or RMS values/shape/WCS are invalid.
    TypeError
        If RMS is an unsupported map representation.
    astropy.units.UnitsError
        If RMS units cannot convert to the cube data unit.

    Notes
    -----
    Errors count finite included channels in the selected slab and assume
    independent channels with constant noise along each sightline. Empty
    sightlines have NaN moment-zero errors. Signal-to-noise ratios are NaN
    where denominators are nonfinite or nonpositive. Inputs and their
    headers are not modified and no derived attributes are attached.
    The helpers materialize data in memory. Without velocity selection,
    any supported spectral unit is retained in moments and errors.

    For included finite intensities I_i at coordinates v_i, let A=sum(I_i),
    mu=mom1, s=sigma, and r=RMS. First-order errors are
    r/abs(A) * sqrt(sum((v_i-mu)**2)) and
    r/(2*s*abs(A)) * sqrt(sum(((v_i-mu)**2-s**2)**2)). These account for
    the shared normalization and estimated centroid. They condition on a
    fixed mask and exclude channel covariance, mask-selection uncertainty,
    calibration errors, and uncertainty in RMS. Zero/nonfinite A or RMS
    and nonfinite required moments give NaN; zero dispersion gives NaN
    sigmaerr. A single-channel centroid has zero error, so its ratio is NaN.
    mom1s2n depends on the coordinate zero point and is not emission-detection
    significance. No positive-intensity selection is imposed.
    """
    if mom_velocity is None or (isinstance(mom_velocity, str) and mom_velocity == ''):
        slab = cube
    else:
        try:
            bounds = u.Quantity(mom_velocity)
            if bounds.unit == u.dimensionless_unscaled:
                bounds = bounds.value * u.km/u.s
            bounds = bounds.to(u.km/u.s)
        except (TypeError, ValueError, u.UnitsError) as exc:
            raise ValueError('Velocity bounds must be two finite velocities.') from exc
        if bounds.shape != (2,) or not np.all(np.isfinite(bounds.value)) or bounds[0] > bounds[1]:
            raise ValueError('Velocity bounds must be two finite, ordered velocities.')
        if not cube.spectral_axis.unit.is_equivalent(u.m/u.s):
            raise ValueError('Convert the cube spectral axis to velocity before selecting a velocity slab.')
        axis = cube.spectral_axis.to(bounds.unit)
        if bounds[1] < axis.min() or bounds[0] > axis.max():
            raise ValueError('Velocity range does not intersect the cube.')
        slab = cube.spectral_slab(*bounds)
    mom0 = slab.moment(order=0, axis=0)
    maximum = slab.max(axis=0)
    error = _moment0_error(slab, rms, mom0)
    rms_values = _rms_quantity(rms, slab.unit, slab.shape[1:], slab.wcs.celestial)
    rms_header = maximum.hdu.header.copy()
    rms_header['BUNIT'] = slab.unit.to_string(format='fits')
    rms_hdu = fits.PrimaryHDU(rms_values.value, rms_header)
    out = dict.fromkeys(['mom0', 'max', 'mom1', 'sigma', 'mom0err', 'maxerr',
                         'mom1err', 'sigmaerr', 'mom0s2n', 'maxs2n', 'mom1s2n', 'sigmas2n'])
    out.update(mom0=mom0.hdu, max=maximum.hdu, mom0err=error.hdu, maxerr=rms_hdu)
    ratios = [('mom0s2n', mom0.value, error.value),
              ('maxs2n', maximum.value, rms_values.value)]
    if velomoms:
        mom1 = slab.moment(order=1, axis=0)
        sigma = slab.linewidth_sigma()
        mom1err, sigmaerr = _velocity_errors(slab, rms_values, mom1, sigma)
        out.update(mom1=mom1.hdu, sigma=sigma.hdu,
                   mom1err=mom1err, sigmaerr=sigmaerr)
        ratios.extend([('mom1s2n', mom1.value, mom1err.data),
                       ('sigmas2n', sigma.value, sigmaerr.data)])
    for name, numerator, denominator in ratios:
        ratio = np.full(mom0.shape, np.nan)
        np.divide(numerator, denominator, out=ratio,
                  where=np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0))
        header = mom0.hdu.header.copy()
        header['BUNIT'] = ''
        out[name] = fits.PrimaryHDU(ratio, header)
    return out
