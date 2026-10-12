"""SpectralCube extensions for noise, moments, and expanded masks."""
from numbers import Real

import numpy as np
import astropy.units as u
from radio_beam import Beam
from spectral_cube import SpectralCube, DaskSpectralCube, VaryingResolutionSpectralCube

from .spectralcubeplus_rms import get_rms_auto as _get_rms_auto, get_rms as _get_rms, _rms_product
from .spectralcubeplus_masks import get_expmask as _get_expmask, _beam_geometry
from .spectralcubeplus_moments import _moment0_error, _rms_quantity

# Register once, rather than changing the unit registry on every read.
_IRAM_UNITS = [u.def_unit('K (Tmb)', represents=u.K),
               u.def_unit('K (Ta*)', represents=u.K)]
u.add_enabled_units(_IRAM_UNITS)


def _print_rms(label, rms):
    """Report a spatial RMS map without printing the whole array."""
    finite = rms.value[np.isfinite(rms.value)]
    if finite.size:
        print('{}: median={:.6g}, range=[{:.6g}, {:.6g}] {}.'.format(
            label, np.median(finite), finite.min(), finite.max(), rms.unit))
    else:
        print('{}: no finite values ({}).'.format(label, rms.unit))


class SpectralCubePlus(SpectralCube):
    """A SpectralCube with spatial noise maps and emission-mask workflows.

    Parameters
    ----------
    *args, **kwargs
        Forwarded to ``spectral_cube.SpectralCube``. Data arrays use
        (spectral, y, x) order; supply a WCS and a data unit as for the parent.

    Attributes
    ----------
    rms : astropy.units.Quantity or None
        Estimated noise map of shape (y, x), in the cube data unit.
    rms_hdu : astropy.io.fits.PrimaryHDU or None
        The same noise values with celestial WCS, BUNIT, and available beam.
    masked : SpectralCubePlus or None
        Cube produced by ``get_expmask``, with RMS products attached.
    mask_cube : SpectralCubePlus or None
        Latest detection cube, with detection RMS attached. Its beam is
        enlarged when ``get_expmask(smooth_factor=...)`` is used.
    expmask : numpy.ndarray
        Boolean (spectral, y, x) inclusion mask, created by ``get_expmask``.
    mom0 : spectral_cube.Projection or None
        Most recent ``moment0(axis=0)`` result.
    mom0err : spectral_cube.Projection or None
        Most recent ``moment0err()`` result. This attribute is distinct from
        the callable method ``moment0err``.

    Notes
    -----
    Cube transformations create new instances and reset derived state. Compute
    noise after preparing the cube, or supply an aligned RMS explicitly to error
    helpers. ``get_expmask`` attaches RMS to its new masked cube explicitly.
    Only common-beam cubes are supported by the reader. Reading with Dask keeps
    that backend, but the added RMS, mask, and error helpers materialize data.
    """
    def __init__(self, *args, **kwargs):
        """Initialize the parent cube and empty derived-product attributes.

        Parameters
        ----------
        *args, **kwargs
            Forwarded unchanged to the parent cube constructor.
        """
        super().__init__(*args, **kwargs)
        self.rms = None
        self.rms_hdu = None
        self.masked = None
        self.mask_cube = None
        self.mom0 = None
        self.mom0err = None

    @classmethod
    def read(cls, filename, *args, **kwargs):
        """Read a cube while preserving its header, mask, beam, and backend.

        Parameters
        ----------
        filename : str, path-like, or file-like
            Input accepted by ``SpectralCube.read``, normally a FITS cube.
        *args, **kwargs
            Forwarded to the parent reader, including ``format`` and ``use_dask``.

        Returns
        -------
        cube : SpectralCubePlus
            Cube with empty derived state. With ``use_dask=True``, the returned
            subclass also inherits from ``DaskSpectralCube``.

        Raises
        ------
        NotImplementedError
            If the input has per-channel beams. Convolve to a common beam first.

        Notes
        -----
        The module registers IRAM ``K (Tmb)`` and ``K (Ta*)`` as Kelvin aliases.
        Reader metadata, spectral-unit scaling, and Dask scheduler settings are
        preserved; errors from the upstream reader propagate to the caller.
        """
        obj = super().read(filename, *args, **kwargs)
        if isinstance(obj, VaryingResolutionSpectralCube) or hasattr(obj, 'beams'):
            raise NotImplementedError('Per-channel beam cubes are not supported; convolve to a common beam first.')
        result_class = _DaskSpectralCubePlus if isinstance(obj, DaskSpectralCube) else cls
        cube = result_class(data=obj._data, wcs=obj.wcs, mask=obj.mask, meta=obj.meta.copy(),
                            header=obj._header.copy(), fill_value=obj._fill_value,
                            beam=obj._meta.get('beam'),
                            allow_huge_operations=obj.allow_huge_operations,
                            wcs_tolerance=obj._wcs_tolerance)
        cube._spectral_unit = obj._spectral_unit
        cube._spectral_scale = obj._spectral_scale
        if isinstance(obj, DaskSpectralCube):
            cube._scheduler_kwargs = obj._scheduler_kwargs.copy()
        return cube

    def get_rms(self, rms_velocities, mask=None):
        """Estimate and store a noise map from selected spectral intervals.

        Parameters
        ----------
        rms_velocities : array-like or astropy.units.Quantity
            One ordered pair, shape (2,), or multiple pairs, shape (N, 2).
            Numeric bounds use the current spectral-axis unit; Quantity bounds
            must convert to that unit. Bounds are inclusive and overlaps count
            each channel once. Choose intervals representative of line-free noise.
        mask : spectral_cube.MaskBase or boolean array, optional
            Additional inclusion mask accepted by ``with_mask``. True includes a
            sample. It is combined with the existing mask for noise estimation
            without changing this cube's mask.

        Returns
        -------
        rms : astropy.units.Quantity
            Gaussian-scaled median absolute deviation (MAD) over selected finite,
            included samples, shape (y, x), in the cube data unit. Empty
            sightlines have NaN noise values.

        Raises
        ------
        ValueError
            If intervals are malformed, nonfinite, reversed, incompatible with
            the spectral unit, or select no channel centres.

        Notes
        -----
        Stores ``rms`` and ``rms_hdu`` and clears the cached ``mom0err`` product.
        Existing expanded masks and masked cubes are not recomputed: regenerate
        these products if you change the noise estimate.
        """
        cube = self if mask is None else self.with_mask(mask)
        self.rms, self.rms_hdu = _get_rms(cube, rms_velocities)
        self.mom0err = None
        return self.rms

    def get_rms_auto(self, mask=None):
        """Estimate and store noise using negative samples and clipping.

        Parameters
        ----------
        mask : spectral_cube.MaskBase or boolean array, optional
            Additional inclusion mask accepted by ``with_mask``, combined with
            the existing cube mask only for this calculation.

        Returns
        -------
        rms : astropy.units.Quantity
            Noise map of shape (y, x), in the cube data unit. All-invalid
            sightlines return NaN.

        Notes
        -----
        Starts with finite-sample, Gaussian-scaled MAD. Sightlines with at least
        10 negative samples are refined using reflected negative samples and a
        channel-count-dependent upper clipping threshold. This assumes noise is
        approximately symmetric about zero and negative samples trace noise;
        absorption or baseline offsets can bias the estimate. See the standalone
        ``get_rms_auto`` for the prescription. Data are loaded into memory.
        Stores ``rms`` and ``rms_hdu`` and clears ``mom0err``. Existing expanded
        masks and masked cubes are not recomputed.
        """
        cube = self if mask is None else self.with_mask(mask)
        self.rms, self.rms_hdu = _get_rms_auto(cube)
        self.mom0err = None
        return self.rms

    def get_expmask(self, *args, rms=None, smooth_factor=None, **kwargs):
        """Create an expanded emission mask and an RMS-equipped masked cube.

        Parameters
        ----------
        rms : Quantity, Projection, or FITS image HDU, optional
            Explicit native-resolution RMS map. Otherwise reuse ``self.rms``
            or estimate it automatically if missing. Mutually exclusive with
            ``smooth_factor``. RMS sources and finite summaries are printed.
        smooth_factor : float, optional
            Finite factor >= 1 multiplying both original beam FWHM axes;
            position angle is preserved. Build the mask using a spatially
            convolved cube and newly estimated automatic detection RMS.
            Factor 1 skips convolution but still estimates detection RMS.
        *args, **kwargs
            Forwarded after the detection cube and RMS to the standalone
            ``spectralcubeplus_masks.get_expmask``. The supported options are:

            hthresh : float, default 5
                Strict high threshold in RMS multiples for emission seeds.
            lthresh : float, default 2
                Strict low threshold in RMS multiples for seed growth.
                Require finite ``0 <= lthresh <= hthresh``.
            beamarea : float, default 2
                Minimum connected area per channel in Gaussian beam areas.
                Finite and nonnegative; zero disables area removal.
            radfrac : float, default 1
                Closing radius as a fraction of half the beam major FWHM.
                Finite and nonnegative; zero gives a one-voxel closing kernel.
            npix : int, default 1
                Seed must also be included at offsets +/-npix channels.
                Nonnegative; zero disables this pruning. No boundary wrapping.

        Returns
        -------
        expmask : numpy.ndarray of bool
            Inclusion array of shape (spectral, y, x), also stored as
            ``self.expmask``. True means included, not masked out.

        Raises
        ------
        ValueError
            If RMS is misaligned, the beam or spatial WCS is invalid, mask
            settings or the smoothing factor are invalid, or explicit RMS
            and a smoothing factor are both supplied.
        TypeError
            If the stored RMS has an unsupported representation.
        astropy.units.UnitsError
            If RMS units cannot convert to the cube data unit.

        Notes
        -----
        Stores native RMS in ``self.rms``/``self.rms_hdu`` and attaches it to
        ``self.masked`` for full-resolution moment uncertainties. Detection
        RMS is attached separately to ``self.mask_cube``. Smoothing preserves
        pixel and spectral grids, so no mask regridding is needed. Pruning and
        closing use the detection beam. ``self.masked`` always contains the
        original-resolution data and beam. This cube's data/mask are unchanged.

        High-threshold seeds are spectrally pruned, grown through connected
        low-threshold emission, closed with a dynamic 3-D kernel, then spatially
        pruned per channel. Closing may include samples below the low threshold.
        The final mask excludes existing masked samples, invalid data, and
        nonfinite or nonpositive RMS. See the standalone helper for geometry and
        connectivity details. This calculation loads the cube into memory.
        """
        if rms is not None and smooth_factor is not None:
            raise ValueError('rms and smooth_factor cannot both be supplied.')
        target_beam = None
        if smooth_factor is not None:
            if (isinstance(smooth_factor, (bool, np.bool_))
                    or not isinstance(smooth_factor, Real)
                    or not np.isfinite(smooth_factor) or smooth_factor < 1):
                raise ValueError('smooth_factor must be a finite number >= 1.')
            _beam_geometry(self, 0)
            target_beam = Beam(major=self.beam.major * smooth_factor,
                               minor=self.beam.minor * smooth_factor,
                               pa=self.beam.pa)

        native_rms = self.rms if rms is None else rms
        native_hdu = self.rms_hdu if rms is None else None
        if native_rms is None:
            print('No stored RMS; calculating native RMS automatically with get_rms_auto().')
            native_rms, native_hdu = _get_rms_auto(self)
            source = 'Automatic native RMS'
        else:
            source = 'Using stored native RMS' if rms is None else 'Using explicit native RMS'
        normalized = _rms_quantity(native_rms, self.unit, self.shape[1:], self.wcs.celestial)
        if rms is not None or not isinstance(native_rms, u.Quantity):
            native_rms = normalized
        if native_hdu is None:
            _, native_hdu = _rms_product(self, normalized.value)
        _print_rms(source, normalized)

        valid = np.isfinite(self.unmasked_data[:].value)
        if self.mask is not None:
            valid &= self.mask.include()
        # Excluded samples must not become finite signal during convolution,
        # even when the original cube uses a finite fill value.
        detection = self.with_mask(valid).with_fill_value(np.nan)
        if target_beam is not None:
            print('Retaining native RMS for uncertainties; determining a new smoothed RMS '
                  'automatically with get_rms_auto().')
            if smooth_factor > 1:
                detection = detection.convolve_to(target_beam)
            detection.get_rms_auto()
            _print_rms('Smoothed detection RMS', detection.rms)
        else:
            detection.rms = native_rms
            detection.rms_hdu = native_hdu.copy()
        expmask = _get_expmask(detection, detection.rms, *args, **kwargs) & valid
        masked = self.with_mask(expmask)

        self.rms = native_rms
        self.rms_hdu = native_hdu
        self.mom0err = None
        self.mask_cube = detection
        self.expmask = expmask
        self.masked = masked
        self.masked.expmask = self.expmask
        self.masked.rms = self.rms
        self.masked.rms_hdu = self.rms_hdu.copy()
        return self.expmask

    def moment0(self, axis=0, how='auto', **kwargs):
        """Integrate along an axis using the inherited moment-zero method.

        Parameters
        ----------
        axis : int, default 0
            Integration axis in (spectral, y, x) order. Axis 0 gives the usual
            spectral integrated-intensity map.
        how : str, default 'auto'
            Upstream computation strategy, e.g. 'auto', 'cube', 'slice', or 'ray'.
        **kwargs
            Additional arguments forwarded to ``SpectralCube.moment0``.

        Returns
        -------
        moment : spectral_cube.Projection
            Integrated map with the data unit times the integrated-axis unit.
            For axis 0, it is also cached as ``self.mom0``.

        Notes
        -----
        Respects the current cube mask. Calling this on the original cube does
        not automatically use ``self.masked``; use ``self.masked.moment0()`` for
        expanded-mask integration. Non-spectral results are not cached in mom0.
        """
        result = super().moment0(axis=axis, how=how, **kwargs)
        if axis == 0:
            self.mom0 = result
        return result

    def moment0err(self, rms=None):
        """Calculate spectral moment-zero uncertainty for this cube's mask.

        Parameters
        ----------
        rms : astropy.units.Quantity, spectral_cube.Projection, or FITS image HDU, optional
            Spatial noise map of shape (y, x), convertible to the cube data unit.
            Defaults to ``self.rms``. A FITS HDU must declare BUNIT. Available
            celestial WCS must match the cube; plain Quantities have no WCS check.

        Returns
        -------
        error : spectral_cube.Projection
            Map of ``RMS * abs(channel_width) * sqrt(N)``, where N counts finite,
            included samples per sightline. Units match spectral moment zero.
            Empty sightlines and nonfinite or zero RMS give NaN. Also stored in
            ``self.mom0err`` without replacing the callable method.

        Raises
        ------
        ValueError
            If RMS is missing, negative, misaligned, or spectral spacing is
            nonuniform or has invalid width.
        TypeError
            If RMS is not a supported map type.
        astropy.units.UnitsError
            If RMS units are incompatible with the cube.

        Notes
        -----
        Assumes independent channels and constant noise along each sightline.
        Does not require an expanded mask or a previously calculated moment.
        Supply RMS explicitly after cube transformations, which reset derived
        state. Counts and moments are calculated from the actual cube mask;
        ``self.masked`` is not selected automatically. Loads data into memory.
        """
        self.mom0err = _moment0_error(self, self.rms if rms is None else rms)
        return self.mom0err

    def writeto(self, filename, *args, **kwargs):
        """Write the cube using the inherited writer.

        Parameters
        ----------
        filename : str or path-like
            Destination accepted by ``SpectralCube.write``.
        *args, **kwargs
            Forwarded to ``write``, including ``overwrite`` and ``format``.

        Returns
        -------
        result
            The upstream writer's return value (normally None for FITS).

        Notes
        -----
        Writes the cube represented by this object, not its attached RMS or
        moment maps. Write those products separately using their HDU or writer.
        """
        return self.write(filename, *args, **kwargs)


class _DaskSpectralCubePlus(SpectralCubePlus, DaskSpectralCube):
    """Combine SpectralCubePlus products with upstream Dask cube operations.

    Notes
    -----
    Selected internally by ``SpectralCubePlus.read(use_dask=True)``. Keeps
    Dask data and scheduler settings; added noise, mask, and channel-count
    helpers still materialize data in memory.
    """
