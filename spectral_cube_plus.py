from spectral_cube import SpectralCube
from spectral_cube import Projection
import astropy.units as u
import numpy as np

from spectralcubeplus_rms import get_rms_auto as _get_rms_auto, get_rms as _get_rms
from spectralcubeplus_masks import get_expmask as _get_expmask
from spectralcubeplus_moments import get_mom0err, get_nchanbool

class SpectralCubePlus(SpectralCube):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # Initialize parent class
        self.rms = None
        self.rms_hdu = None
        self.masked = None
        self.mom0 = None
        self.mom0err = None
        self.masked = None
    
    @classmethod
    def read(cls, filename, *args, **kwargs):

        # Add the IRAM units to the enabled units
        u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
        u.add_enabled_units(u.def_unit(['K (Ta*)'], represents=u.K))
        
        # Use the parent class's read method to create an instance
        obj = super().read(filename, *args, **kwargs)
        # Create an instance of SpectralCubePlus using the SpectralCube object
        return cls(data=obj._data, wcs=obj._wcs, mask=obj._mask, meta=obj._meta)

    def get_rms(self, mask=None, *args, **kwargs):
        if mask is not None:
            self.rms, self.rms_hdu = _get_rms(self.with_mask(mask), *args, **kwargs)
        else:
            self.rms, self.rms_hdu = _get_rms(self, *args, **kwargs)  # Use the inherited cube object
        return self.rms

    def get_rms_auto(self, mask=None):
        if mask is not None:
            self.rms, self.rms_hdu = _get_rms_auto(self.with_mask(mask))
        else:
            self.rms, self.rms_hdu = _get_rms_auto(self)  # Use the inherited cube object
        return self.rms
    
    def get_expmask(self, *args, **kwargs): 
        expmask = _get_expmask(self, self.rms, *args, **kwargs)

        def _set_mask(mask):
            return mask == True
        
        self.expmask = expmask
        self.masked = self.with_mask(expmask)
        
        self.masked.expmask = expmask
        self.masked.rms = self.rms
        self.masked.rms_hdu = self.rms_hdu

        return self.expmask

    def moment0(self):
        self.mom0 = self.moment(order=0)
        return self.mom0

    def moment0err(self):
        
        nchan = get_nchanbool(self.expmask)
        cdelt3 = np.abs(self.hdu.header['CDELT3']) * self.spectral_axis.unit
        mom0 = self.mom0

        if mom0 is None:
            mom0 = self.masked.mom0

        mom0err = get_mom0err(mom0, nchan, self.rms_hdu, cdelt3)
        self.mom0err = Projection(mom0err.data, wcs=mom0.wcs, unit=mom0.unit)

        return self.mom0err

    def writeto(self, filename, *args, **kwargs):
        super().write(filename, *args, **kwargs)