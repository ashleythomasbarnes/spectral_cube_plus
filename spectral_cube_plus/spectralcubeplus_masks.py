"""Threshold expansion and spatial/spectral pruning."""
import numpy as np
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_area, proj_plane_pixel_scales
from scipy import ndimage
from spectral_cube.utils import NoBeamError
from .spectralcubeplus_moments import _rms_quantity


def get_threshmask(cube, thresh_map, thresh=0):
    """Compare a cube with a scaled threshold map.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Cube in (spectral, y, x) order.
    thresh_map : array-like or astropy.units.Quantity
        Threshold values, normally shape (y, x), broadcastable against the
        cube. Numeric values use the cube data unit; Quantities must have
        compatible units. This helper does not check map WCS alignment.
    thresh : float, default 0
        Multiplier applied to the threshold map. The comparison is strict
        ``cube > thresh_map * thresh``.

    Returns
    -------
    mask : spectral_cube.MaskBase
        Comparison mask, evaluated as a boolean inclusion array via
        ``mask.include()``.
    masked_cube : spectral_cube.SpectralCube
        New cube with the comparison applied in addition to its existing
        mask. The input cube is not modified.
    """
    if not isinstance(thresh_map, u.Quantity):
        thresh_map = np.asarray(thresh_map) * cube.unit
    mask = cube > thresh_map * thresh
    return mask, cube.with_mask(mask)


def get_circmask(h=11, w=11, l=11, center=None, radius=5):
    """Create a spherical boolean structuring element in array coordinates.

    Parameters
    ----------
    h, w, l : int, default 11
        Lengths of axes 0, 1, and 2, respectively. Coordinates are voxel
        indices, not angular or spectral world coordinates.
    center : float or array-like of length 3, optional
        Sphere centre. A scalar is used on all axes; a three-vector gives
        coordinates in axis order. None or the legacy empty string uses
        ``(shape - 1) / 2`` separately on each axis.
    radius : float, default 5
        Radius in voxel units. Samples exactly on the sphere are included.

    Returns
    -------
    structure : numpy.ndarray of bool
        Array of shape (h, w, l), True inside the sphere. A sphere extending
        beyond these dimensions is clipped to the array.

    Notes
    -----
    The expanded-mask helper builds its own beam-scaled kernel rather than
    using this fixed-size element.
    """
    if center is None or (isinstance(center, str) and center == ''):
        center = (np.array([h, w, l]) - 1) / 2
    center = np.broadcast_to(center, (3,))
    indices = np.ogrid[:h, :w, :l]
    distance2 = sum((coordinate - origin)**2 for coordinate, origin in zip(indices, center))
    return distance2 <= radius**2


def get_prunemask(mask, thresh):
    """Remove small spatial components independently in each channel.

    Parameters
    ----------
    mask : array-like
        Three-dimensional inclusion mask in (spectral, y, x) order; values
        are converted to bool.
    thresh : float
        Finite, nonnegative minimum component size in spatial pixels.
        Components with size equal to the threshold are retained; zero
        disables removal of included components.

    Returns
    -------
    pruned : numpy.ndarray of bool
        Pruned copy with the input shape. The input is not modified.

    Raises
    ------
    ValueError
        If the threshold is invalid or the mask is not three-dimensional.

    Notes
    -----
    Uses four-connected components in each 2-D spatial plane. Diagonally
    adjacent pixels alone are not connected. There is no component linking
    between spectral channels.
    """
    if not np.isfinite(thresh) or thresh < 0:
        raise ValueError('Spatial pruning threshold must be finite and nonnegative.')
    result = np.asarray(mask, dtype=bool).copy()
    if result.ndim != 3:
        raise ValueError('Mask must be three-dimensional.')
    for k, plane in enumerate(result):
        labels, _ = ndimage.label(plane)
        sizes = np.bincount(labels.ravel())
        keep = sizes >= thresh
        keep[0] = False
        result[k] = keep[labels]
    return result


def get_prunemaskvelo(mask, npix=1):
    """Keep seeds with included samples at a channel offset on both sides.

    Parameters
    ----------
    mask : array-like
        Three-dimensional inclusion mask in (spectral, y, x) order.
    npix : int, default 1
        Nonnegative channel offset. A sample survives only if it and the
        samples at +/-npix are included. Zero returns a boolean copy without
        pruning; bool values are not accepted as integers.

    Returns
    -------
    pruned : numpy.ndarray of bool
        New array with the input shape. With npix > 0, the first and last
        npix channels are excluded. If no channel has both neighbours,
        the result is entirely False. The input is unchanged.

    Raises
    ------
    ValueError
        If npix is not a nonnegative integer or mask is not three-dimensional.

    Notes
    -----
    Does not wrap at boundaries or require every intervening channel when
    npix > 1. This is an offset check, not a minimum consecutive-run length.
    """
    if isinstance(npix, (bool, np.bool_)) or not isinstance(npix, (int, np.integer)) or npix < 0:
        raise ValueError('npix must be a nonnegative integer.')
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 3:
        raise ValueError('Mask must be three-dimensional.')
    if npix == 0:
        return mask.copy()
    result = np.zeros_like(mask)
    if 2 * npix < mask.shape[0]:
        result[npix:-npix] = mask[npix:-npix] & mask[:-2*npix] & mask[2*npix:]
    return result


def _beam_geometry(cube, radfrac):
    """Derive Gaussian beam area and a beam-scaled 3-D closing kernel.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Common-beam cube with finite positive celestial pixel scales/area.
    radfrac : float
        Nonnegative radius multiplier, validated by the caller.

    Returns
    -------
    beam_pixels : float
        ``pi * major * minor / (4 * ln(2) * pixel_area)``, in spatial pixels.
    structure : numpy.ndarray of bool
        Odd-sized ellipsoidal kernel in (spectral, y, x) order. Spatial radii
        are half the beam major FWHM times radfrac, divided by each pixel
        scale. The channel radius uses the smallest spatial pixel scale,
        not a spectral resolution. Zero radius returns a one-voxel kernel.

    Raises
    ------
    ValueError
        If the beam or celestial pixel geometry is missing or invalid.
    """
    try:
        major = cube.beam.major.to_value(u.deg)
        minor = cube.beam.minor.to_value(u.deg)
    except (AttributeError, ValueError, NoBeamError) as exc:
        raise ValueError('Expanded masking requires a valid spatial beam.') from exc
    spatial_wcs = cube.wcs.celestial
    scales = np.abs(proj_plane_pixel_scales(spatial_wcs))
    area = abs(proj_plane_pixel_area(spatial_wcs))
    if (len(scales) != 2 or not np.all(np.isfinite(scales)) or np.any(scales <= 0)
            or not np.isfinite(area) or area <= 0 or not np.isfinite(major)
            or not np.isfinite(minor) or major <= 0 or minor <= 0):
        raise ValueError('Beam and spatial pixel scales must be finite and positive.')
    beam_pixels = np.pi * major * minor / (4 * np.log(2) * area)
    # Retain a 3-D closing kernel; x/y radii represent the same angular size.
    # Its channel radius follows the major beam radius in the smallest pixels.
    angular_radius = major / 2 * radfrac
    radii = np.array([angular_radius / scales.min(),
                      angular_radius / scales[1], angular_radius / scales[0]])
    if angular_radius == 0:
        return beam_pixels, np.ones((1, 1, 1), dtype=bool)
    half_sizes = np.ceil(radii).astype(int)
    grid = np.ogrid[tuple(slice(-n, n+1) for n in half_sizes)]
    structure = sum((axis / radius)**2 for axis, radius in zip(grid, radii)) <= 1
    return beam_pixels, structure


def get_expmask(cube, rms, hthresh=5, lthresh=2, beamarea=2, radfrac=1, npix=1):
    """Grow high-threshold emission seeds, close gaps, and prune small areas.

    Parameters
    ----------
    cube : spectral_cube.SpectralCube
        Input cube in (spectral, y, x) order, with a valid common spatial
        beam and celestial WCS. Its existing mask is respected.
    rms : astropy.units.Quantity, spectral_cube.Projection, or FITS image HDU
        Noise map of shape (y, x), convertible to the cube data unit.
        FITS PrimaryHDU/ImageHDU inputs must declare BUNIT. Available
        celestial WCS must match the cube. Plain arrays and scalars are not
        accepted. Nonfinite and zero RMS sightlines are excluded.
    hthresh : float, default 5
        Strict high threshold in RMS multiples for initial seeds.
    lthresh : float, default 2
        Strict low threshold in RMS multiples for allowed growth. Require
        finite ``0 <= lthresh <= hthresh``.
    beamarea : float, default 2
        Minimum connected area per channel, in Gaussian beam areas. Finite
        and nonnegative. Zero disables spatial area removal.
    radfrac : float, default 1
        Closing radius relative to half the beam major FWHM. Finite and
        nonnegative; zero disables closing via a one-voxel kernel.
    npix : int, default 1
        High-threshold seeds must have high-threshold samples at offsets
        +/-npix channels. Nonnegative; zero disables this pruning. No wrapping
        or requirement on intervening channels for npix > 1.

    Returns
    -------
    mask : numpy.ndarray of bool
        Inclusion array with the cube shape. True means included. Neither
        the cube nor the RMS input is modified; no products are attached.

    Raises
    ------
    ValueError
        If RMS is missing, negative, or misaligned; beam/pixel geometry is
        invalid; or thresholds, radius, area, or channel offset are invalid.
    TypeError
        If RMS is not a Quantity, Projection, or FITS image HDU.
    astropy.units.UnitsError
        If RMS units are incompatible with the cube.

    Notes
    -----
    Thresholding uses finite, currently included samples and strictly
    positive finite RMS. Spectrally pruned high seeds grow to convergence
    through the low mask using six-connected 3-D dilation. Binary closing
    then uses a dynamically sized ellipsoidal kernel. Its angular spatial
    radius is ``beam.major / 2 * radfrac``; its channel radius is that angular
    radius divided by the smallest spatial pixel scale, not a velocity width.
    Default closing treats samples outside the cube as False and can remove
    edge emission. Closing may add samples below lthresh. The result is
    intersected with valid samples before four-connected spatial components
    smaller than ``beamarea * pi * major * minor / (4 * ln(2) * pixel_area)``
    are removed separately in each channel. Spectral seed pruning is not
    reapplied after growth/closing. The full cube is loaded into memory.
    """
    if rms is None:
        raise ValueError('Compute RMS first or supply an RMS map.')
    if (not np.all(np.isfinite([hthresh, lthresh, beamarea, radfrac]))
            or lthresh < 0 or hthresh < lthresh or beamarea < 0 or radfrac < 0):
        raise ValueError('Require 0 <= lthresh <= hthresh and nonnegative beamarea/radfrac.')
    rms = _rms_quantity(rms, cube.unit, cube.shape[1:], cube.wcs.celestial)
    beam_pixels, structure = _beam_geometry(cube, radfrac)
    data = cube.unmasked_data[:].value
    valid = np.isfinite(data)
    if cube.mask is not None:
        valid &= cube.mask.include()
    valid &= np.isfinite(rms.value)[None, :, :] & (rms.value > 0)[None, :, :]
    low = valid & (data > lthresh * rms.value)
    high = valid & (data > hthresh * rms.value)
    high = get_prunemaskvelo(high, npix=npix)
    grown = ndimage.binary_dilation(high, iterations=-1, mask=low)
    closed = ndimage.binary_closing(grown, structure=structure) & valid
    return get_prunemask(closed, thresh=beam_pixels * beamarea)
