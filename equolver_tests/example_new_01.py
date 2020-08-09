#!/usr/bin/env python
import warnings
import copy
import numpy as np
from scipy import stats
import radio_beam
from astropy.io import fits
from astropy import constants
from astropy import wcs
from astropy import units as u
from kapteyn import wcs as kwcs

class Beach:
    """
    Lame acronym for a utility to equalize synthesized beams over
    channels

    Class variables:
    HIFREQ: rest frequency of the HI line
    
    """
    HIFREQ = 1420405751.7667

    def __init__(self, cubenames = None,
                       bmaj = np.nan, bmaj_replace = False,
                       bmin = np.nan, bmin_replace = False,
                       bpa = np.nan, bpa_replace = False,
                       restfreq = HIFREQ, restfreq_replace = False,
                       normfreq = 1E9,
                       parameter = 'all',
                       stype = ['median','madstdev'], sample = 'all', percents = 90,
                       tolerance = 0.1, nsamps = 200, epsilon = 0.0005,
                       tar_bmaj_inter = ['bmaj', 'median', 'total'],
                       tar_bmaj_slope = ['bmaj', 'madstdev', 'total'],
                       tar_bmaj_absc = 3.0,
                       tar_bmin_inter = ['bmaj', 'median', 'total'],
                       tar_bmin_slope = ['bmaj', 'madstdev', 'total'],
                       tar_bmin_absc = 3.0,
                       tar_bpa_inter =  ['bpa', 'median', 'total'],
                       tar_bpa_slope =  ['bpa', 'madstdev', 'total'],
                       tar_bpa_absc =  0.0,
                       tar_scaling = 'freq', genbstats_exe = True, gentarget_exe = True,
                       gentrans_exe = True, transnames = None, overwrite = True,
                       verb = False):
        """
        Private instance variables:
        (multiple: None, a float, a list of floats, a numpy array, 
         or a list with numpy arrays)

        _cubenames (list)     : Names of input cubes
        _headers (list)         : List of open header hdus

        _bmaj     (multiple)  : Beam major axis default value(s)
        _bmaj_replace (bool)  : Enforce usage of default values?
                                  (True = yes)
        _bmin     (multiple)  : Beam minor axis default value(s)
        _bmin_replace (bool)  : Enforce usage of default values?
                              (True = yes)
        _bpa     (multiple)   : Beam position angle default value(s)
        _bpa_replace (bool)   : Enforce usage of default values? 

        _restfreq (multiple)  : Rest frequency default value(s)
        _bpa_replace (bool)   : Enforce usage of default values? 

        _binfo_input (list)   : List of arrays of point spread func-
                                tion, in the order of headernames and
                                headers slowest index same as index of
                                header in headers, followed by _bmaj,
                                _bmin, and _bpa, fastest index
                                channels, so it is an array of size
                                (len(_headernames),3,chans), where
                                chans is the number of channels.
        _binfo_pixel (list)  : _binfo_input converted into pixel
                               scaling using dispersion instead of HPBW
        _bstats (dict)       : Dictionary containing all statistics
        _binfo_target (list)  : Target beam properties
        """
        self._initvars()
        self._verb = verb
        self.initheaders(cubenames = cubenames)
        self._initbinfo_inputvar(bmaj = bmaj, bmaj_replace =
                                 bmaj_replace, bmin = bmin,
                                 bmin_replace = bmin_replace,
                                 bpa = bpa, bpa_replace =
                                 bpa_replace, restfreq =
                                 restfreq, restfreq_replace =
                                 restfreq_replace, verb =
                                 self._verb)
        
        if genbstats_exe:
            self._genbstats_exe = True

        self.genbinfo(verb = self._verb)

        self._normfreq = copy.deepcopy(normfreq)

        for para in ['parameter', 'stype', 'sample',
                     'percents', 'tolerance', 'nsamps', 'epsilon',
                     'gentarget_exe']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])

        if self._genbstats_exe:
            self.genbstats(verb = self._verb)

        for para in [ 'tar_bmaj_inter', 'tar_bmaj_slope',
                      'tar_bmaj_absc', 'tar_bmin_inter',
                      'tar_bmin_slope', 'tar_bmin_absc',
                      'tar_bpa_inter', 'tar_bpa_slope', 'tar_bpa_absc',
                      'tar_scaling','gentrans_exe']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])
            
        if self._gentarget_exe:
            self.gentarget(verb = self._verb)
            
        for para in [ 'gentrans_exe', 'transnames', 'overwrite']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])
            
        if self._gentrans_exe:
            self.gentrans(verb = self._verb)
            
        return

    def _initvars(self):
        """
        Reset/init instance variables
        """

        # Headers
        self._cubenames = None
        self._headers = None

        # Defaults
        self._bmaj = None
        self._bmaj_replace = False
        self._bmin = None
        self._bmin_replace = False
        self._bpa = None
        self._bpa_replace = False
        self._restfreq = None
        self._restfreq_replace = False

        # This is the normalisation for frequency-dependent beam
        # properties
        self._normfreq = None

        # bmaj, bmin, pa, nu, pixelsize
        self._binfo_input = None

        # bmaj, bmin, sin pa, cos pa, nu
        self._binfo_pixel = None

        self._parameter = None
        self._scaling = None
        self._stype = None
        self._sample = None
        self._percents = None
        self._tolerance = None
        self._nsamps = None
        self._epsilon = None

        # Statistics
        self._bstats = None
        
        self._tar_bmaj_inter = None
        self._tar_bmaj_slope = None
        self._tar_bmaj_absc = None
        self._tar_bmin_inter = None
        self._tar_bmin_slope = None
        self._tar_bmin_absc = None
        self._tar_bpa_inter = None
        self._tar_bpa_slope = None
        self._tar_bpa_absc = None
        self._tar_scaling = None

        # bmaj, bmin, sin pa, cos pa, nu
        self._binfo_target = None

        self._transnames = None
        self._overwrite = None
        
        return

    @property
    def cubenames(self):
        """
        Return a copy of cubenames
        """
        if type(self._cubenames) == type(None):
            return None
        return copy.deepcopy(self._binfo_input)

    @cubenames.setter
    def cubenames(self, value):

        # There is no use in letting the user change the cube names
        # but not changing the cubes, so this is enforced
        self.initheaders(cubenames = value)
        return
    
    @cubenames.deleter
    def cubenames(self, value):
        self._resetheaders()
        
    @property
    def headers(self):
        """
        Return a copy of the list headers

        The headers themselves are not copies
        """
        # Notice that this is not a deep copy to save memory
        # But this is really the responsibility of the user
        
        if type(self._headers) == type(None):
            return None
        return copy.copy(self._headers)

    def _returndefault(self, value):
        if type(value) == type(None):
            return None
        returnlist = []

        # This is a list of np arrays
        for item in value:
            returnlist.append(item.copy())
        return returnlist
    
    @property
    def bmaj(self):
        """
        Return a copy of bmaj
        """
        return copy.deepcopy(self._bmaj)

    @bmaj.setter
    def bmaj(self, value):
        """
        Set bmaj
        """
        self._bmaj = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @bmaj.deleter
    def bmaj(self, value):
        self._bmaj = np.nan
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def bmaj_replace(self):
        return self._bmaj_replace

    @bmaj_replace.setter
    def bmaj_replace(self, value):
        self._bmaj_replace = value
    
    @bmaj_replace.deleter
    def bmaj_replace(self, value):
        self._bmaj_replace = False

    @property
    def bmin(self):
        """
        Return a copy of bmin
        """
        return copy.deepcopy(self._bmin)

    @bmin.setter
    def bmin(self, value):
        """
        Set bmin
        """
        self._bmin = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @bmin.deleter
    def bmin(self, value):
        self._bmin = np.nan
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def bmin_replace(self):
        return self._bmin_replace

    @bmin_replace.setter
    def bmin_replace(self, value):
        self._bmin_replace = value
    
    @bmin_replace.deleter
    def bmin_replace(self, value):
        self._bmin_replace = False

    @property
    def bpa(self):
        """
        Return a copy of bpa
        """
        return copy.deepcopy(self._bpa)

    @bpa.setter
    def bpa(self, value):
        """
        Set bpa
        """
        self._bpa = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @bpa.deleter
    def bpa(self, value):
        self._bpa = np.nan
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def bpa_replace(self):
        return self._bpa_replace

    @bpa_replace.setter
    def bpa_replace(self, value):
        self._bpa_replace = value
    
    @bpa_replace.deleter
    def bpa_replace(self, value):
        self._bpa_replace = False

    @property
    def restfreq(self):
        """
        Return a copy of restfreq
        """
        return copy.deepcopy(self._restfreq)

    @restfreq.setter
    def restfreq(self, value):
        """
        Set restfreq
        """
        self._restfreq = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @restfreq.deleter
    def restfreq(self, value):
        self._restfreq = None
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def restfreq_replace(self):
        return self._restfreq_replace

    @restfreq_replace.setter
    def restfreq_replace(self, value):
        self._restfreq_replace = value
    
    @restfreq_replace.deleter
    def restfreq_replace(self, value):
        self._restfreq_replace = False

    @property
    def normfreq(self):
        """
        Return a copy of normfreq
        """
        return copy.deepcopy(self._normfreq)

    @normfreq.setter
    def normfreq(self, value):
        """
        Set normfreq
        """
        self._normfreq = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @normfreq.deleter
    def normfreq(self, value):
        self._normfreq = 1.0E9
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def parameter(self):
        """
        Return a copy of parameter
        """
        return copy.deepcopy(self._parameter)

    @parameter.setter
    def parameter(self, value):
        """
        Set parameter
        """
        self._parameter = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @parameter.deleter
    def parameter(self, value):
        self._parameter = 'all'
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def stype(self):
        """
        Return a copy of stype
        """
        return copy.deepcopy(self._stype)

    @stype.setter
    def stype(self, value):
        """
        Set stype
        """
        self._stype = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @stype.deleter
    def stype(self, value):
        self._stype = ['median','medstdev']
        self._bstats = None
        self._binfo_target = None
        return
    
    @property
    def sample(self):
        """
        Return a copy of sample
        """
        return copy.deepcopy(self._sample)

    @sample.setter
    def sample(self, value):
        """
        Set sample
        """
        self._sample = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @sample.deleter
    def sample(self, value):
        self._sample = 'all'
        self._bstats = None
        self._binfo_target = None
        return
    
    @property
    def percents(self):
        """
        Return a copy of percents
        """
        return copy.deepcopy(self._percents)

    @percents.setter
    def percents(self, value):
        """
        Set percents
        """
        self._percents = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @percents.deleter
    def percents(self, value):
        self._percents = 90
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def tolerance(self):
        """
        Return a copy of tolerance
        """
        return copy.deepcopy(self._tolerance)

    @tolerance.setter
    def tolerance(self, value):
        """
        Set tolerance
        """
        self._tolerance = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @tolerance.deleter
    def tolerance(self, value):
        self._tolerance = 0.01
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def nsamps(self):
        """
        Return a copy of nsamps
        """
        return copy.deepcopy(self._nsamps)

    @nsamps.setter
    def nsamps(self, value):
        """
        Set nsamps
        """
        self._nsamps = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @nsamps.deleter
    def nsamps(self, value):
        self._nsamps = 200
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def epsilon(self):
        """
        Return a copy of epsilon
        """
        return copy.deepcopy(self._epsilon)

    @epsilon.setter
    def epsilon(self, value):
        """
        Set epsilon
        """
        self._epsilon = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @epsilon.deleter
    def epsilon(self, value):
        self._epsilon = 0.0005
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def tar_bmaj_inter(self):
        """
        Return a copy of tar_bmaj_inter
        """
        return copy.deepcopy(self._tar_bmaj_inter)

    @tar_bmaj_inter.setter
    def tar_bmaj_inter(self, value):
        """
        Set tar_bmaj_inter
        """
        self._tar_bmaj_inter = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bmaj_inter.deleter
    def tar_bmaj_inter(self, value):
        self._tar_bmaj_inter = ['bmaj', 'freq', 'median', 'all']
        self._binfo_target = None
        return

    @property
    def tar_bmaj_slope(self):
        """
        Return a copy of tar_bmaj_slope
        """
        return copy.deepcopy(self._tar_bmaj_slope)

    @tar_bmaj_slope.setter
    def tar_bmaj_slope(self, value):
        """
        Set tar_bmaj_slope
        """
        self._tar_bmaj_slope = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bmaj_slope.deleter
    def tar_bmaj_slope(self, value):
        self._tar_bmaj_slope = ['bmaj', 'freq', 'madstdev', 'all']
        self._binfo_target = None
        return

    @property
    def tar_bmaj_absc(self):
        """
        Return a copy of tar_bmaj_absc
        """
        return copy.deepcopy(self._tar_bmaj_absc)

    @tar_bmaj_absc.setter
    def tar_bmaj_absc(self, value):
        """
        Set tar_bmaj_absc
        """
        self._tar_bmaj_absc = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bmaj_absc.deleter
    def tar_bmaj_absc(self, value):
        self._tar_bmaj_absc = 3.0
        self._binfo_target = None
        return

    
    @property
    def tar_bmin_inter(self):
        """
        Return a copy of tar_bmin_inter
        """
        return copy.deepcopy(self._tar_bmin_inter)

    @tar_bmin_inter.setter
    def tar_bmin_inter(self, value):
        """
        Set tar_bmin_inter
        """
        self._tar_bmin_inter = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bmin_inter.deleter
    def tar_bmin_inter(self, value):
        self._tar_bmin_inter = ['bmaj', 'freq', 'median', 'all']
        self._binfo_target = None
        return

    @property
    def tar_bmin_slope(self):
        """
        Return a copy of tar_bmin_slope
        """
        return copy.deepcopy(self._tar_bmin_slope)

    @tar_bmin_slope.setter
    def tar_bmin_slope(self, value):
        """
        Set tar_bmin_slope
        """
        self._tar_bmin_slope = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bmin_slope.deleter
    def tar_bmin_slope(self, value):
        self._tar_bmin_slope = ['bmaj', 'freq', 'madstdev', 'all']
        self._binfo_target = None
        return

    @property
    def tar_bmin_absc(self):
        """
        Return a copy of tar_bmin_absc
        """
        return copy.deepcopy(self._tar_bmin_absc)

    @tar_bmin_absc.setter
    def tar_bmin_absc(self, value):
        """
        Set tar_bmin_absc
        """
        self._tar_bmin_absc = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bmin_absc.deleter
    def tar_bmin_absc(self, value):
        self._tar_bmin_absc = 3.0
        self._binfo_target = None
        return

    @property
    def tar_bpa_inter(self):
        """
        Return a copy of tar_bpa_inter
        """
        return copy.deepcopy(self._tar_bpa_inter)

    @tar_bpa_inter.setter
    def tar_bpa_inter(self, value):
        """
        Set tar_bpa_inter
        """
        self._tar_bpa_inter = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bpa_inter.deleter
    def tar_bpa_inter(self, value):
        self._tar_bpa_inter = ['bpa', 'freq', 'median', 'all']
        self._binfo_target = None
        return

    @property
    def tar_bpa_slope(self):
        """
        Return a copy of tar_bpa_slope
        """
        return copy.deepcopy(self._tar_bpa_slope)

    @tar_bpa_slope.setter
    def tar_bpa_slope(self, value):
        """
        Set tar_bpa_slope
        """
        self._tar_bpa_slope = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bpa_slope.deleter
    def tar_bpa_slope(self, value):
        self._tar_bpa_slope = ['bpa', 'freq', 'madstdev', 'all']
        self._binfo_target = None
        return

    @property
    def tar_bpa_absc(self):
        """
        Return a copy of tar_bpa_absc
        """
        return copy.deepcopy(self._tar_bpa_absc)

    @tar_bpa_absc.setter
    def tar_bpa_absc(self, value):
        """
        Set tar_bpa_absc
        """
        self._tar_bpa_absc = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_bpa_absc.deleter
    def tar_bpa_absc(self, value):
        self._tar_bpa_absc = 3.0
        self._binfo_target = None
        return

    @property
    def tar_scaling(self):
        """
        Return a copy of tar_scaling
        """
        return copy.deepcopy(self._tar_scaling)

    @tar_scaling.setter
    def tar_scaling(self, value):
        """
        Set tar_scaling
        """
        self._tar_scaling = copy.deepcopy(value)
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @tar_scaling.deleter
    def tar_scaling(self, value):
        self._tar_scaling = 'freq'
        self._binfo_target = None
        return

    @property
    def genbstats_exe(self):
        """
        Return a copy of genbstats_exe
        """
        return self._genbstats_exe

    @genbstats_exe.setter
    def genbstats_exe(self, value):
        """
        Set genbstats_exe
        """
        self._genbstats_exe = value
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @genbstats_exe.deleter
    def genbstats_exe(self):
        self._genbstats_exe = False
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def gentarget_exe(self):
        """
        Return a copy of gentarget_exe
        """
        return self._gentarget_exe

    @gentarget_exe.setter
    def gentarget_exe(self, value):
        """
        Set gentarget_exe
        """
        self._gentarget_exe = value
        if self._gentarget_exe:
            self.gentarget(verb = False)
        return

    @gentarget_exe.deleter
    def gentarget_exe(self):
        self._gentarget_exe = False
        self._binfo_target = None
        return
    
    @property
    def binfo_input(self):
        """
        Return a copy of binfo_input
        """
        return self._returndefault(self._binfo_input)

    @binfo_input.setter
    def binfo_input(self, value):
        """
        Copy the input
        """
        if type(value) == type(None):
            self._binfo_input = None
            return

        # Input must be a list of np arrays, which we will copy
        returnlist = []
        for i in value:
            returnlist.append(i.copy())
        self._binfo_input = returnlist
        return

    @binfo_input.deleter
    def binfo_input(self, value):
        self._binfo_input = None
    
    @property
    def binfo_pixel(self):
        """
        Return a copy of binfo_pixel
        """
        return self._returndefault(self._binfo_pixel)

    @property
    def binfo_target(self):
        """
        Return a copy of binfo_target
        """
        return self._returndefault(self._binfo_target)

    @property
    def bstats(self):
        """
        Return a copy of bstats
        """
        return copy.deepcopy(self._bstats)

    @property
    def verb(self):
        """
        Return a copy of verb
        """
        return self._verb

    @verb.setter
    def verb(self, value):
        """
        Set verb
        """
        self._verb = value
        return

    @verb.deleter
    def verb(self):
        self._verb = False
        return

    def _unitsconv(self, quantity, units = None):
        """
        Convert to intrinsic units (deg and Hz)

        If unit is specified, try only that unit
        """
        if units == None:
            units = [u.deg, u.Hz]
        
        if type(quantity) == type(1.*u.deg):
            searching = True
            for i in units:
                try:
                    outquant = quantity.to(i).value
                    searching = False
                    break
                except u.core.UnitConversionError:
                    pass
            if searching:
                raise(u.core.UnitConversionError('Incompatible unit'))
        else:
            outquant = quantity
            
        return outquant

    def _resetheaders(self):
        """
        Close all headers in instance and set headers and headernames to None
        """
        if type(self._headers) != type(None):
            self._headers = None
        self._cubenames = None
        return

    def initheaders(self, cubenames = None, silent = False):
        """
        Add cubenames to intrinsic cubenames

        Input:
        cubenames (list of str): List of names of input cubes/images

        Reads cubenames in as the list of target headers/images
        """
        if type(cubenames) == type(None):
            return

        # Check if cubenames have been defined before
        #dontinit = True
        #for i in range(len(cubenames)):
        #    if cubenames[i] != self._cubenames[i]:
        #        dontinit == False
        #        break
            
        # Save the time to reload the cubes    
        #if dontinit and silent:
        #    return
        
        self._resetheaders()
        self._cubenames = copy.deepcopy(cubenames)
        self._headers = []

        if type(self._cubenames) == type(''):
            cubenamelist = [self._cubenames]
        else:
            cubenamelist = self._cubenames
            
        for cube in self._cubenames:
            opencube = fits.open(cube)
            self._headers += [opencube[0].header]
            opencube.close()

        # Cascade down
        self.genbinfo(verb = self._verb)
        return
    
    def _initbinfo_input(self):
        """
        Initiate the beam info arrays for the cubes and fill with nans

        Will initiate a list of ndarrays, one for each header in 
        self.headers with the dimension nchan x 5, where nchan is the
        number of channels in the corresponding cube. Each value is
        set to np.nan
        """
        
        if type(self._headers) == type(None):
            self.initheaders()
            
        if type(self._headers) == type(None):
            return

        self._binfo_input = []
        
        for header in self._headers:
            
            # Create array of nans, size naxis3 x 5 (bmaj, bmin, bpa,
            # psize, freq) and append to the info list
            if 'NAXIS3' in header.keys():
                naxis3 = header['NAXIS3']
            else:
                naxis3 = 1
        

            self._binfo_input.append(np.empty((naxis3,5,)))
            self._binfo_input[-1][:] = np.nan
        return

    def _getdefault(self, inquant, quantname = None):
        """
        Fill target quantname, wich is the name of an instance
        variable, with the content of inquant or return formatted
        inquant
        
        Input:
        (multiple: None, a float, a list of floats, a numpy array, 
         or a list with numpy arrays)
        inquant (multiple): input
        quantname (str)   : Name of quantity (without leading _)

        self._quantname is assumed to be a list of linear ndarrays of
        the size of naxis3 of self._headers (i.e. the number of
        channels), which should be filled with the default values for
        bmaj, bmin, etc., following the same expansion scheme. If
        inquant is None, the call is ignored, if it is a single
        (float) number (including np.nan), all fields in the output
        are filled with that number, if it is a list of numbers
        (including np.nan) each ndarray in the output list is filled
        with the corresponding number, if it is an ndarray, every
        array in the output list is filled with this ndarray, and
        finally, if a list of ndarrays is provided, the output will
        completely be filled with the input. An input list has always
        to have as many elements as self._headers and an ndarray always
        has to have the same dimension as the naxis3 of the
        corresponding cube.

        """

        cinquant = copy.deepcopy(inquant)
        
        if type(cinquant) == type(None) and \
           quantname != None and \
           type(self.__dict__['_'+quantname]) == type(None):
            cinquant = np.nan
        
        if type(cinquant) == type(None):
            return

        # Invoke initbinfo_input
        if type(self._binfo_input) == type(None):
            self._initbinfo_input()

        # If binfo_input is empty, return
        if type(self._binfo_input) == type(None):
            return

        output = []

        for i in range(len(self._headers)):
            
            cinquant = self._unitsconv(cinquant)
            
            # One number for all
            if type(cinquant) == type(.1) or type(cinquant) == type(np.nan):
                output.append(np.repeat(
                    cinquant, self._binfo_input[i].shape[0]))

            # If input is an ndarray it is a default array with
            # one cinquant per channel for all cubes
            elif type(cinquant) == type(np.array([])): 
                output.append(np.repeat(
                    cinquant, self._binfo_input[i].shape[0]))

            # If it is a list, there are two possibilities
            elif type(cinquant) == type([]):

                cinquant = self._unitsconv(cinquant)
                
                # One number for all channels
                if type(cinquant[i]) == type(.1) or \
                type(cinquant[i]) == type(np.nan):
                    output.append(np.repeat(
                        cinquant[i], self._binfo_input[i].shape[0]))

                # The only other possibility is an ndarray
                elif type(cinquant[i]) == type(np.array([])):
                    output.append(cinquant[i])

                # Else doom
                else:
                    raise BaseException('{:s} must be None, a float, a list of',
                                        'floats, a numpy array, or a list with',
                                        'numpy arrays.'.format(quantname))

            # Else doom
            else:
                raise BaseException('{:s} must be None, a float, a numpy',
                                    'array, or a list with numpy arrays.'
                                    .format(quantname))

        if type(quantname) != type(None):
            self.__dict__['_'+quantname] = output
        
        return output

    def _setdefreplace(self, bmaj_replace = None,
                            bmin_replace = None,
                            bpa_replace = None,
                            restfreq_replace = None):
        """
        Read in values for replacements of defvalues
        """
        
        if self._bmaj_replace == None:
            self._bmaj_replace = False
        if self._bmin_replace == None:
            self._bmin_replace = False
        if self._bpa_replace == None:
            self._bpa_replace = False
        if self._restfreq_replace == None:
            self._restfreq_replace = False

        if bmaj_replace != None:
            self._bmaj_replace = bmaj_replace
        if bmin_replace != None:
            self._bmin_replace = bmin_replace
        if bpa_replace != None:
            self._bpa_replace = bpa_replace
        if restfreq_replace != None:
            self._restfreq_replace = restfreq_replace
            
    def _initbinfo_inputvar(self, bmaj = None, bmaj_replace = None, bmin = None,
                    bmin_replace = None, bpa = None, bpa_replace = None,
                            restfreq = None, restfreq_replace = None, verb = False):
        """
        Check existence of variables, return True if a parameter is ill defined
        """
        output = False
        paras = locals().copy()
        paras.pop('self')
        for param in paras.keys():
            if type(paras[param]) != type(None):
                self.__dict__['_'+param] = copy.deepcopy(paras[param])
            else:
                if self.__dict__['_'+param] == None:
                    if self._verb or verb:
                        warnings.warn('Parameter {} is not defined.'.format(param))
                    output = True
        return output
    
    def genbinfo(self, bmaj = None, bmaj_replace = None,
                       bmin = None, bmin_replace = None,
                       bpa = None,  bpa_replace = None,
                       restfreq = HIFREQ, restfreq_replace = None,
                       verb = True):
        """
        Fill beam properties into info table
        
        Input:
        (multiple: None, a float, a list of floats, a numpy
                 array, or a list with numpy arrays)

        bmaj     (multiple): Beam major axis default value(s)
        bmaj_replace (bool): Enforce usage of default values?
                             (True = yes)
        bmin     (multiple): Beam minor axis default value(s)
        bmin_replace (bool): Enforce usage of default values?
                             (True = yes)
        bpa     (multiple) : Beam position angle default value(s)
        bpa_replace (bool) : Enforce usage of default values? 

        restfreq (multiple): Rest frequency default value(s)
        bpa_replace (bool) : Enforce usage of default values? 

        Note that if None is passed as a value, the input is ignored
        if self._quantity is not None. I that case
        """
        stop = self._initbinfo_inputvar(bmaj = bmaj,
                                        bmaj_replace = bmaj_replace,
                                        bmin = bmin, bmin_replace =
                                        bmin_replace, bpa = bpa,
                                        bpa_replace = bpa_replace,
                                        restfreq = restfreq,
                                        restfreq_replace =
                                        restfreq_replace, verb = verb)

        if stop:
            if verb or self._verb:
                warnings.warn('Parameters missing. Not generating beam info.')
            return

        if type(self._headers) == type(None):
            if verb or self._verb:
                warnings.warn('No headers loaded. Returning.')
            return
        
        if type(self._binfo_input) == type(None):
            self._initbinfo_input()
            
        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warning('binfo_input not present. Use initheaders first.')
            return
        
        # Get defaults
        bmaj = self._getdefault(self._bmaj)
        bmin = self._getdefault(self._bmin)
        bpa = self._getdefault(self._bpa)
        restfreq = self._getdefault(self._restfreq)

        #self._setdefreplace(bmaj_replace = bmaj_replace,
                            #bmin_replace = bmin_replace,
                            #bpa_replace = bpa_replace,
                            #restfreq_replace = restfreq_replace)

        for i in range(len(self._headers)):

            # It would be rather surprising if there is a channel-dependent
            # rest frequency, but who knows
            restfreqtab = self._getchanval(
                'RESTFREQ', self._headers[i], value = restfreq[i],
                usevalue = restfreq_replace, usedefault = True)
            
            # Count numbers of axes
            naxis = self._headers[i]['NAXIS']

            pixcrd = np.column_stack((
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64),
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64)))
            pixcrd2 = np.column_stack((
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64),
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64)+1.))

            stocol = 3
            if naxis > 2:
                if self._headers[i]['CTYPE3'] != 'STOKES':
                    pixcrd = np.column_stack((
                        pixcrd, np.arange(self._binfo_input[i].shape[0],
                                          dtype=np.float64)+1.))
                    pixcrd2 = np.column_stack((
                        pixcrd2, np.arange(self._binfo_input[i].shape[0],
                                           dtype=np.float64)+1.))
                else:
                    stocol = 2
            else:
                stocol = 2
                
            for j in range(stocol,naxis):             
                pixcrd = np.column_stack((
                    pixcrd, np.ones(self._binfo_input[i].shape[0],
                                    dtype=np.float64)))
                pixcrd2 = np.column_stack((
                    pixcrd2, np.ones(self._binfo_input[i].shape[0],
                                     dtype=np.float64)))

            # Use Kapteyn to find out about frequency and pixel size
            worldcrd = []
            worldcrd2 = []
            wcshand2 = None
            for j in range(self._binfo_input[i].shape[0]):
                if not np.isnan(restfreqtab[j]):
                    self._headers[i]['RESTFREQ'] = restfreqtab[j]
                wcshand = kwcs.Projection(self._headers[i])
                if stocol == 3:
                    wcshand2 = wcshand.spectra('FREQ')
                    worldcrd += [wcshand2.toworld(pixcrd[j])[2]]
                    worldcrd2 += [np.fabs(
                    wcshand2.toworld(pixcrd[j])[1]-
                    wcshand2.toworld(pixcrd2[j])[1])]
                else:
                    worldcrd2 += [np.fabs(
                    wcshand.toworld(pixcrd[j])[1]-
                    wcshand.toworld(pixcrd2[j])[1])]

            if worldcrd != []:
                self._binfo_input[i][:,3] = np.array(worldcrd)
            self._binfo_input[i][:,4] = np.array(worldcrd2)

            # Now check for cellscal and change the pixel size
            # accordingly
#            print(repr(self._cubes[i][0].header))
            if 'CELLSCAL' in self._headers[i].keys():
                if self._headers[i]['CELLSCAL'] == '1/F':
                    if not type(wcshand2) == type(None):
                        if not np.isnan(self._binfo_input[i][:,3].sum()):
                            self._binfo_input[i][:,4] = wcshand2.crval[2]* \
                                self._binfo_input[i][:,4]/ \
                                self._binfo_input[i][:,3]

            cellscal_use_constant = True
            if type(wcshand2) == type(None):
                cellscal_use_constant = False
            if 'CELLSCAL' in self._headers[i].keys():
                if self._headers[i]['CELLSCAL'] == '1/F':
                    cellscal_use_constant = False
            
            # Determine the scaling of the beam properties
            # Spectral cube and constant cells means that the beam
            # changes reciprocal to frequency
            dscal = 1.
            if cellscal_use_constant:
                if not type(wcshand2) == type(None):
                    if not np.isnan(self._binfo_input[i][:,3].sum()):
                        dscal = wcshand2.crval[2]/ \
                            self._binfo_input[i][:,3]

            # finlist: slowest index cube in 'headers' list, followed
            # by bmaj, bmin, and bpa
            self._binfo_input[i][:,0] = self._getchanval( \
                'BMAJ', self._headers[i], value = bmaj[i],
                usevalue = bmaj_replace, usedefault = True, dscal = dscal)
            self._binfo_input[i][:,1] = self._getchanval( \
                'BMIN', self._headers[i], value = bmin[i],
                usevalue = bmin_replace, usedefault = True, dscal = dscal)
            self._binfo_input[i][:,2] = self._getchanval( \
                'BPA', self._headers[i], value = bpa[i],
                usevalue = bpa_replace, usedefault = True)

        self._genbinfo_pixel(verb = verb)

        if self._genbstats_exe:
            self.genbstats(verb = verb)
        return

    def _getchanval(self, prefix, thedict, value = None, usevalue = False,
                   usedefault = True, dscal = 1.):
        """
        Return an ndarray of values from an input dict with formatted
        keywords

        Input:
        prefix (string)         : Prefix of formatted keywords
        thedict (dict)          : Input dict
        value (ndarray)         : Default values
        usevalue (None or bool) : If set output is identical to
                                   value (or np.nan if value is None)
        dscal (float or ndarray): Scale default with this (channel-
                                   dependent) value
        
    
        Rather specialized method. thedict is assumed to be a
        dictionary-like representing a fits header. Then a linear
        ndarray with length 'NAXIS3' is created, initially filled with
        nans. If 'NAXIS3' does not exist, the length of the array is
        1. If usevalue is True, then the output array is filled with
        the values as specified for the parameter value (which) is
        assumed to be a linear ndarray of the same length as the
        output.  If not, the method assumes that there is a number of
        keywords with the same prefix prefix, e.g. 'prefix1',
        'prefix2', ..., in the dict, where the number after the prefix
        stands for a corresponding channel. Then the contents of the
        output array are assigned the corres- ponding values for the
        given channels (first channel to be 1 in the key as per FITS
        standard, hence having an index 0 in the array).  If
        usedefault is True, the input dict will be searched for a key
        'prefix' without extension, which is then used as a default
        should an extended keyword not exist for a channel. value is
        used as a further default should that key not be found in the
        input dict (independently of usedefault; set value = None in
        order not to use it as a default). If a default value is used,
        it gets multiplied by dscale, which is either a scalar or an
        ndarray. If it is the latter, each value gets multiplied indi-
        vidually.

        Example: 
        ourdict = { 'BMAJ1': 0.1, 'BMAJ3': 0.2, 'BMAJ': 0.3, 'NAXIS3': 3}
        _getchanval('BMAJ', ourdict, value = 17.3, usevalue = False, 
                   usedefault = True)
        [0.1, 0.3, 0.2]

        ourdict = { 'BMAJ1': 0.1, 'BMAJ3': 0.2, 'NAXIS3': 3}
        _getchanval('BMAJ', ourdict, value = 17.3, usevalue = False,
                   usedefault = True)
        [0.1, 17.3, 0.2]

        _getchanval('BMAJ', ourdict, value = [17.2,17.3,17.4], 
                   usevalue = True, usedefault = True)
        [17.2, 17.3, 17.4]

        _getchanval('BMAJ', ourdict, value = None, usevalue = False, 
                   usedefault = False)
        [0.1, nan, 0.2]

        """
        # Check if there is actually a third dimension
        if 'NAXIS3' in thedict.keys():
            naxis3 = thedict['NAXIS3']
        else:
            naxis3 = 1
        
        # Start with an empty array
        output = np.empty((naxis3,))
        output[:] = np.nan
        
        if usevalue:
            if value != None:
                output[:] = value
                    
            return output
            
        # Go through dict and search for a cubedefault
        cubedefault = np.nan
        if usedefault:
            for key in thedict.keys():
                if key == prefix:
                    cubedefault = thedict[prefix]
                    break

        # 'Channel' numbers
        chnum = []

        # Values
        chval = []
    
        # Go through dict and assign numbers
        for key in thedict.keys():
            if key.startswith(prefix):
                if key != prefix:
                    try:
                        chnum = int(key[len(prefix):])
                        try:
                            chval = float(thedict[key])
                            output[chnum-1] = chval
                        except:
                            pass
                    except:
                        pass

        # Now fill remaining nans
        output[np.isnan(output)] = cubedefault

        # Scale output
        output[np.isnan(output)] = (output*dscal)[np.isnan(output)]
        
        if type(value) != type(None):
            if type(value) == type(np.array([])):
                output[np.isnan(output)] = value[np.isnan(output)]
            else:
                output[np.isnan(output)] = value
        
        return output

    def _genbinfo_input(self, verb = False):
        """
        Alias for use in genbstats
        """
        self.genbinfo(verb = verb)
        return
    
    def _genbinfo_pixel(self, verb = False):
        """
        Generate binfo_pixel, a list of ndarrays, which is
        the beam info file binfo_input converted to pixel
        units.
        """

        
        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warnings.warn('Trying to generate beam info.')
            self.genbinfo()
            
        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warnings.warn('No beam information read, which disables'+ \
                              ' further processing.')
            return
        
        for binfarray in self._binfo_input:
            if binfarray[:,0].sum() == np.nan:
                if verb or self._verb:
                    warnings.warn('Not sufficient information about beam major'+ \
                                  ' axes, which disables further processing.')
                return
        
            if binfarray[:,1].sum() == np.nan:
                if verb or self._verb:
                    warnings.warn('Not sufficient information about beam minor'+ \
                                  ' axes, which disables further processing.')
                return
            
            if binfarray[:,2].sum() == np.nan:
                if verb or self._verb:
                    warnings.warn('Not sufficient information about beam '+ \
                                  'position angle, which disables further '+ \
                                  'processing.')
                return
            
            if binfarray[:,4].sum() == np.nan:
                if verb or self._verb:
                    warnings.warn('Not sufficient information about pixel '+ \
                                  'dimension, which disables further '+ \
                                  'processing.')
                return
            
        self._binfo_pixel = []
        for binfarray in self._binfo_input:
            boutfarray = binfarray.copy()
            
            # Scale to pixels as unit and dispersion
            boutfarray[:,0] = binfarray[:,0]/binfarray[:,4]/np.sqrt(np.log(256))
            boutfarray[:,1] = binfarray[:,1]/binfarray[:,4]/np.sqrt(np.log(256))

            # Scale to cos pa and sin pa
            boutfarray[:,2] = np.pi*binfarray[:,2]/180.
            boutfarray[:,3] = np.sin(boutfarray[:,2])
            boutfarray[:,4] = np.cos(boutfarray[:,2])

            # Copy frequency
            # boutfarray[:,4] = binfarray[:,3]
            
            self._binfo_pixel.append(boutfarray)
        return
        
    def _initbstats(self, verb = False):
        """
        Create a structure able to contain all statistics

        The intrinsic statistics structure _bstats contains all
        statistics. It is a nested dictionary. See genbstats for a
        description. This method generates the structure and fills it
        with nans.

        """

        # Count the number of cubes and the maximum number of channels
        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warnings.warn('Trying to generate beam statistics.')
            self.genbinfo()

        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warnings.warn('No beam information read in, which disables'+ \
                              ' further processing.')
            return

        ncubes = len(self._binfo_input)

        channels = 0
        for i in range(ncubes):
            channels = max(self._binfo_input[i].shape[0], channels)

        if channels == 0:
            if verb or self._verb:
                warnings.warn('No cube with more than 0 channels present. '+ \
                              ' This is virtually impossible. No statistics '+ \
                              'structure built.')
            return
            
        self._bstats = {}
        # parameter
        for key0 in ['bmaj', 'bmin', 'bpa']:
            self._bstats[key0] = {}
            
            # scaling
            #for key1 in ['input', 'freq']:
            #    self._bstats[key0][key1] = {}

            # stype
            for key1 in ['minimum', 'maximum', 'average', 'stdev', 'median',
                         'mad', 'madstdev', 'percentile', 'percents',
                         'commonbeam']:
                self._bstats[key0][key1] = {}
                # for key4 in ['cube', 'chan', 'tota']:

                # sample
                self._bstats[key0][key1]['total'] = np.nan
                self._bstats[key0][key1]['chan'] = \
                    np.empty((channels,))
                self._bstats[key0][key1]['chan'][:] = np.nan
                self._bstats[key0][key1]['cube'] = []
                for i in range(ncubes):
                    self._bstats[key0][key1]['cube'].append(np.nan)
                        
        return
    
    def _initbstatsvar(self, parameter = None, stype = None,
                       sample = None, percents = None, tolerance = None,
                       nsamps = None, epsilon = None, verb = False):
        """
        Check existence of variables, return True if a parameter is ill defined
        """
        output = False
        paras = locals().copy()
        paras.pop('self')
        for param in paras.keys():
            if type(paras[param]) != type(None):
                self.__dict__['_'+param] = copy.deepcopy(paras[param])
            else:
                if self.__dict__['_'+param] == None:
                    if verb or self._verb:
                        warnings.warn('Parameter {} is not defined.'.format(param))
                    output = True
        return output
            
    def genbstats(self, parameter = None, stype = None, sample = None,
                  percents = None, tolerance = None, nsamps = None,
                  epsilon = None, verb = True):
        """
        Generate statistics and dump it into the bstats structure
        
        Input:
        parameter (str or list of str): Parameter name ('bmaj', 
                                        'bmin', 'bpa')
        stype (str or list of str)    : Type of statistics to
                                        calculate ('minimum',
                                        'maximum', 'average',
                                        'stdev', 'median', 'mad',
                                        'madstdev', 'percentile',
                                        'percents', 'combeam')
        sample (str or list of str)   : Sample(s) to calculate
                                        statistics on ('cube', 'chan',
                                        'total')
        perc (float)                  : Percents for the percentile
                                        statistics
        tolerance (float)             : Tolerance for the common beam
        nsamps (int)                  : Number of edges of beam for
                                        common beam
        epsilon (float)               : Epsilon for common beam

        The method generates statistics on the collected beam
        properties. The parameters parameter, scaling, sample, and
        type determine which part of the bstats structure gets
        filled. If for any of the parameters 'all' is chosen (which is
        the default), all fields are filled. The scaling type
        corresponds to the tables/structures binfo_input and
        binfo_pixel generated using methods genbinfo, _genbinfo_pixel.

        From top to bottom level:

        parameter: 
            major axis dispersions/hpbws ('bmaj')
            minor axis dispersions/hpbws ('bmin')
            beam position angles ('bpa')

        stype:
            Minimum ('minimum')
            Maximum ('maximum')
            Average ('average')
            Standard deviation ('stdev')
            Median ('median')
            Median-absolute-deviation ('mad')
            Standard deviation calculated from the
                median-absolute-deviation ('madstdev')
            score at percents ('percentile')
            Percents percents is corresponding to ('percents')

            Common beam as calculated using the radio-beam module
            (https://radio-beam.readthedocs.io) based on the Khachiyan
            algorithm (https://en.wikipedia.org/wiki/Ellipsoid_method).
            Parameters tolerance, nsamps, epsilon are used for this 
            method ('commonbeam')

        sample:
            statistics to be carried out for all channels per cube
            ('cube', generates lists with length of the number of
            input cubes)

            statistics to be carried out for all cubes per channel
            ('chan', generates linear ndarrays with a length of the
            maximum number of channels in any cube)

            statistics to be carried out for all channels in all
            cubes ('total', generates a float)

        """

        stop = self._initbstatsvar(parameter = parameter,
                                stype = stype, sample = sample,
                                percents = percents, tolerance = tolerance,
                                nsamps = nsamps, epsilon = epsilon,
                                verb = verb)

        if stop:
            if verb or self._verb:
                warnings.warn('Parameters missing. Not generating statistics.')
            return

        if self._bstats == None:
            if verb or self._verb:
                warnings.warn('Attempting to initialize bstats')
            self._initbstats()
        
        if self._bstats == None:
            if verb or self._verb:
                warnings.warn('Failing to initialize bstats. Returning.')
            return
        
        para = copy.deepcopy(self._parameter)
        styp = copy.deepcopy(self._stype)
        samp = copy.deepcopy(self._sample)

        # Expand all parameters to the same format
        if type(para) == type(''):
            para = [para]
        if type(styp) == type(''):
            styp = [styp]
        if type(samp) == type(''):
            samp = [samp]

        #para = self._parameter[:]
        #scal = self._scaling[:]
        #styp = self._stype[:]
        #samp = self._sample[:]

        if 'all' in para:
            para = ['bmaj', 'bmin', 'bpa']
        if 'all' in styp:
            styp = ['minimum', 'maximum', 'average', 'stdev', 'median',
                    'mad', 'madstdev', 'percentile', 'percents',
                    'commonbeam']
        if 'all' in samp:
            samp = ['cube', 'chan', 'total']
        
        # Make sure we do have the required information available
        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warnings.warn('Attempting to genereate input'+ \
                              ' information struct.')
            self._genbinfo_input()
        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warnings.warn('Inut information struct '+ \
                              'cannot be generated. Returning without '+ \
                              'generating statistics.')
            return

        # This is needed later, map for parameters
        parma = {'bmaj': 0, 'bmin': 1, 'bpa': 2}
        
        # Now let's do it

        # get struct
        #struct = self._binfo_input

        for sam in samp:

            # collect sample as list of np arrays
            collist = []

            #print(sam)

            if sam == 'cube':
                for i in self._binfo_input:
                    collist.append(i[:,[0,1,2]])
            elif sam == 'chan':
                for i in range(self._bstats['bmaj']['maximum'][sam].shape[0]):
                    collist.append(np.empty((0,3,), dtype=np.float64))
                    for j in self._binfo_input:
                        if j.shape[0] > i:
                            collist[i] = np.append(collist[i],
                                                   j[i:i+1,[0,1,2]],
                                                   axis = 0)
            elif sam == 'total':
                collist = [np.empty((0,3,))]
                for i in self._binfo_input:
                    collist[0] = np.append(collist[0], i[:,[0,1,2]],
                                           axis = 0)
            else:
                if verb or self._verb:
                    warnings.warn('Statistics cannot be generated.'+ \
                                  'Check sample parameter. Returning.')
            #print('collist', collist)

            #stackedcahn = 1
            #if collist[0].shape[0] > stackedcahn:
            #    print(collist[0],collist[0][stackedcahn,2])
            #print('yo\n')
            for sty in styp:
                met = getattr(self, '_gen'+sty)

                # Apply statistics
                kwargs = {
                    'percents' : self._percents,
                    'tolerance': self._tolerance,
                    'nsamps'   : self._nsamps,
                    'epsilon'  : self._epsilon
                    }
                stats = met(collist, **kwargs)
                #print('stats ', stats)

                # Get function to generate statistics
                for par in para:

                    #print('par', par, 'parma', parma[par])

                    # Put the right value into array
                    self._bstats[par][sty][sam] = stats[:,parma[par]]
                    #print('self._bstats[', par, '][', sty, '][', sam,'] =', self._bstats[par][sca][sty][sam])

        # Cascade down
        if self._gentarget_exe:
            self.gentarget(verb = verb)
            
    def _genminimum(self, parray, **kwargs):
        """
        Return an ndarray with minima for bmin, bmaj, and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        kwargs (dict): further arguments

        For each element of the list generates the minimum along the
        first axis for each of the m components of the ndarray.
        kwargs is a dummy dict.
        """
        
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(output, np.amin(array, axis = 0).reshape((1,3,)),
                               axis = 0)

        return output
    
    def _genmaximum(self, parray, **kwargs):
        """
        Return an ndarray with maxima for bmin, bmaj, and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        kwargs (dict): further arguments

        For each element of the list generates the maximum along the
        first axis for each of the m components of the ndarray.
        kwargs is a dummy dict.
        """
        
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(output, np.amax(array, axis = 0).reshape((1,3,)),
                               axis = 0)

        return output
        
    def _genaverage(self, parray, **kwargs):
        """
        Return an ndarray with averages for bmin, bmaj, and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        kwargs (dict): further arguments

        For each element of the list generates the mean along the
        first axis for each of the m components of the ndarray.
        kwargs is a dummy dict.
        """
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(output, np.mean(array, axis = 0).reshape((1,3,)),
                               axis = 0)

        return output
        
    def _genstdev(self, parray, **kwargs):
        """
        Return an ndarray with standard deviations for bmin, bmaj, 
        and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        kwargs (dict): further arguments

        For each element of the list generates the standard deviation
        along the first axis for each of the m components of the
        ndarray.
        kwargs is a dummy dict.
        """

        output = np.empty((0, 3))
        for array in parray:
            
            # Make sure that there is no error if there is only one
            if array.shape[0] < 2:
                dof = 0.
            else:
                dof = 1.
        
            output = np.append(output, np.std(array, axis = 0,
                                      ddof = dof).reshape((1,3,)), axis = 0)

        return output
    
    def _genmedian(self, parray, **kwargs):
        """
        Return an ndarray with medians for bmin, bmaj, and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        kwargs (dict): further arguments

        For each element of the list generates the median along the
        first axis for each of the m components of the ndarray.
        kwargs is a dummy dict.
        """
        
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(output, np.median(array, axis = 0).reshape((1,3,)),
                               axis = 0)

        return output

    def _genmad(self, parray, **kwargs):
        """
        Return an ndarray with median absolute deviations for bmin,
        bmaj, and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        kwargs (dict): further arguments

        For each element of the list generates the median absolute
        deviations along the first axis for each of the m components
        of the ndarray.
        kwargs is a dummy dict.
        """
        
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(output, stats.median_abs_deviation(
                array, axis = 0, scale = 1.).reshape((1,3,)), axis = 0)

        return output

    def _genmadstdev(self, parray, **kwargs):
        """
        Return an ndarray with standard deviation as derived from the
        median absolute deviations for bmin, bmaj, and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)

        For each element of the list generates the standard deviation
        as derived from the median absolute deviations along the first
        axis for each of the m components of the ndarray.
        kwargs is a dummy dict.
        """
        
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(output, stats.median_abs_deviation(
                array, axis = 0, scale = 'normal').reshape((1,3,)), axis = 0)

        return output

    def _genpercentile(self, parray, **kwargs):
        """
        Return an ndarray with scores at percentss for bmin, bmaj,
        and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        perc (float) : Percents used for the calculation

        For each element of the list generates the perc percents
        along the first axis for each of the m components of the ndarray
        """
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(output, np.percentile(
                array, kwargs['percents'], axis = 0).reshape((1,3,)), axis = 0)

        return output

    def _genpercents(self, parray, **kwargs):
        """
        Return an ndarray describing a common beam for bmin, bmaj,
        and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)

        Common beam as calculated using the radio-beam module
        (https://radio-beam.readthedocs.io) based on the Khachiyan
        algorithm (https://en.wikipedia.org/wiki/Ellipsoid_method)
        ('commonbeam')
        kwargs is a dummy argument
        """
        
        output = np.empty((0, 3))
        for array in parray:
            output = np.append(
                output, np.array([kwargs['percents']]).repeat(3).reshape((1,3,)),
                axis = 0)

        return output
        
    def _gencommonbeam(self, parray, **kwargs):
        """
        Return an ndarray with input percents for bmin, bmaj,
        and bpa (or more)
        
        Input:
        parray (list): list of ndarrays of shape (n, m,)
        perc (float) : Percents used for the calculation

        Just dump the input perc to the output
        """
        
        output = np.empty((0, 3))
        for array in parray:

            # Unit of this one does not matter
            bmaj=array[:,0] * u.deg

            # Unit of this one does not matter
            bmin=array[:,1] * u.deg

            # Unit of this one has to be deg
            bpa=array[:,2] * u.deg

            my_beams = radio_beam.Beams(bmaj, bmin, bpa)

            repeat = True
            tolerance = kwargs['tolerance']
            while repeat:
                try:
                    common_beam = my_beams.common_beam(
                        tolerance = tolerance,
                        nsamps = kwargs['nsamps'],
                        epsilon = kwargs['epsilon'])
                    repeat = False
                except radio_beam.utils.BeamError:
                    tolerance = tolerance/1.1
                    
            bmaj = common_beam.major.to(u.degree).value
            bmin = common_beam.minor.to(u.degree).value
            bpa = common_beam.pa.to(u.degree).value
            
            output = np.append(
                output, np.array([bmaj, bmin, bpa]).reshape((1,3,)),
                axis = 0)

        return output
        
    def _getar(self, item):
        """
        Decode item and return part of a beam struct

        Input:
        item (multiple)

        If item is a list of four strings, the function will return a
        binfo struct based on the statistics as saved in bstats, based
        on the three strings, which should list the parameter, scaling,
        stype, and sample that the generated struct should be based
        on. The result is a list with the length of the number of
        input cubes, each element being an ndarray of dimension
        nplanes x 1 , where nplanes is the number of planes per cube.
        Alternatively, the parameter is read in using _getdefault.

        """

        # First check which type of input this is
        if type(item) == type([]):
            if len(item) == 3:
                
                # The alternative is either numeric or quantity,
                # so this is sufficient for identifying the input
                if type(item[0]) == type(''):
                    # print(self._bstats)
                    # Ensure that we can actually do this
                    if type(self._bstats) == type(None) or \
                       self._bstats[item[0]][item[1]][item[2]] \
                       == np.nan:
                        if verb or self._verb:
                            warnings.warn('Attempting to genereate part of '+ \
                                          'beam statistics struct.')
                        self.genbstats(parameter = item[0], stype =
                                       item[1], sample = item[2], verb = verb)
                    if type(self._bstats) == type(None) or \
                       self._bstats[item[0]][item[1]][item[2]] \
                       == np.nan:
                        if verb or self._verb:
                            warnings.warn('Failed to genereate beam statistics'+ \
                                          'struct. Returning empty-handed.')
                        return               

                    # This should have cascaded up, such that the information
                    # is at our hands
                    output = []
                    for i in range(len(self._headers)):
                        output.append(np.zeros(self._binfo_input[i].shape[0]))

                        if item[2] == 'cube':
                            output[i] = output[i]+\
                                self._bstats[item[0]][item[1]][item[2]][i]
                        elif item[2] == 'channel' or item[2] == 'total':
                            output[i] = output[i]+\
                                self._bstats[item[0]][item[1]][item[2]]
                        else:
                            raise BaseException('Argument must be \'cube\', ',
                                                '\'channel\', or \'total\.')
                    return output
                
        # Now we know this can only be numeric and a direct input
        output = self._getdefault(item)
        return output
                
    def _inittargetvar(self, tar_bmaj_inter = None, tar_bmaj_slope = None,
                       tar_bmaj_absc = None, tar_bmin_inter = None,
                       tar_bmin_slope = None, tar_bmin_absc = None,
                       tar_bpa_inter = None, tar_bpa_slope = None,
                       tar_bpa_absc = None, tar_scaling = None,
                       verb = False):
        """
        Check existence of variables, return True if a parameter is ill defined
        """
        output = False
        paras = locals().copy()
        paras.pop('self')
        paras.pop('verb')
        for param in paras.keys():
            if type(paras[param]) != type(None):
                self.__dict__['_'+param] = copy.deepcopy(paras[param])
            else:
                if self.__dict__['_'+param] == None:
                    if verb or self._verb:
                        warnings.warn('Parameter {} is not defined.'.format(param))
                    output = True
        return output

    def gentarget(self, tar_bmaj_inter = None, tar_bmaj_slope = None,
                  tar_bmaj_absc = None, tar_bmin_inter = None,
                  tar_bmin_slope = None, tar_bmin_absc = None,
                  tar_bpa_inter = None, tar_bpa_slope = None,
                  tar_bpa_absc = None, tar_scaling = None,
                  verb = True):
        """Generate a target beam structure

        Input:
        (multiple: None, a float, a list of floats, a numpy array, 
         a list with numpy arrays, or a quadruple of strings)
        tar_bmaj_inter (multiple): Beam major axis intercept
        tar_bmaj_slope (multiple): Beam major axis slope
        tar_bmaj_absc  (multiple): Beam major axis abscissae
        tar_bmin_inter (multiple): Beam minor axis intercept
        tar_bmin_slope (multiple): Beam minor axis slope    
        tar_bmin_absc  (multiple): Beam minor axis abscissae
        tar_bpa_inter  (multiple): Beam position angle axis intercept
        tar_bpa_slope  (multiple): Beam position angle axis slope    
        tar_bpa_absc   (multiple): Beam position angle abscissae
        tar_scaling    (multiple): Use 1/F scaling when calculating
                               the target array, either 'freq' or
                               'input'

        The method first generates a struct binfo_target of the same
        dimension as binfo_input (or, equivalently or
        binfo_pixel). Each element in that struct (a list with ncubes
        ndarrays with dimension channels x 5) represents a target
        quantity bmaj, bmin, bpa, frequency, pixel scale, where the
        last two are copy-pasted from the input beam struct
        (binfo_input), and is generated as follows:
        
        For each quantity bmaj, bmin, bpa the para-
        meters quant_inter (intercept), quant_slope (slope),
        quant_absc (abscissa) result in the generation of a list with
        ndarrays of dimension channels x 1. Then the output column
        for quant in binfo_target gets calculated as: 

        quant = quant_inter + quant_slope*quant_absc
 
        The parameters are either direct inputs of the target
        quantities with multiple possible input as described for
        _getdefault. Alternatively, a list of four strings is passed,
        denoting (in that order) parameter, scaling, stype, and sample
        as described in method genbstats. The corresponding values
        will then be copied from the bstats struct.

        Depending on scaling, the struct is then transformed into an
        analogue of a binfo_pixel struct, which then represents the
        target quantities in pixel coordinates.  Notice that freq
        means the inverse of the transformation of the input beam info
        struct into the frequency-scaled one.  Major and minor axis
        beams entered directly are then interpreted as beams at a
        frequency as specified in the instance variable normfreq
        (which is a list of length of the number of cubes of linear
        ndarrays of length of the single cubes).

        That final struct is stored as member binfo_target of the
        class instance.

        """

        # Initialisation of required arrays will happen in the
        # helper methods, so we do not have to care, apart from this:
        
        stop = self._inittargetvar(tar_bmaj_inter = tar_bmaj_inter,
                                   tar_bmaj_slope = tar_bmaj_slope, tar_bmaj_absc =
                                   tar_bmaj_absc, tar_bmin_inter = tar_bmin_inter,
                                   tar_bmin_slope = tar_bmin_slope, tar_bmin_absc =
                                   tar_bmin_absc, tar_bpa_inter = tar_bpa_inter,
                                   tar_bpa_slope = tar_bpa_slope, tar_bpa_absc =
                                   tar_bpa_absc, tar_scaling = tar_scaling, verb =
                                   verb)
        
        if stop:
            if verb or self._verb:
                warnings.warn('Parameters missing. Not generating target properties.')
            return
                                
        # Do not try to regenerate, if this is none, it has deliberately been
        # set to None, as __init__ has a default of 1E9
        if type(self._normfreq) == type(None):
            if verb or self._verb:
                warnings.warn('No normalization frequency read in, which '+ \
                              'disables further processing.')
            return

        normfreq = self._getdefault(self._normfreq)

        # Check if input is sufficient, checking for bstats is enough
        # as this cascades up
        if type(self._bstats) == type(None):
            if verb or self._verb:
                warnings.warn('Attempting to generate bstats')
            self.genbstats()
        if type(self._bstats) == type(None):
            if verb or self._verb:
                warnings.warn('Failed to generate bstats. Not generating target.')
            return
        tar_bmaj_inter = self._getar(self._tar_bmaj_inter)
        tar_bmaj_slope = self._getar(self._tar_bmaj_slope)
        tar_bmaj_absc  = self._getar(self._tar_bmaj_absc)
        tar_bmin_inter = self._getar(self._tar_bmin_inter)
        tar_bmin_slope = self._getar(self._tar_bmin_slope)
        tar_bmin_absc  = self._getar(self._tar_bmin_absc)
        tar_bpa_inter  = self._getar(self._tar_bpa_inter)
        tar_bpa_slope  = self._getar(self._tar_bpa_slope)
        tar_bpa_absc   = self._getar(self._tar_bpa_absc)

        # For the preliminary output we make a copy of input
        targar = copy.deepcopy(self._binfo_input)

        for i in range(len(targar)):

            # Notice that this is the inverse operation to translating
            # the input beams to frequency scaling
            if tar_scaling == 'freq':
                scale = targar[i][:,3]/normfreq[i]
            else:
                scale = np.ones(targar[i].shape[0])

            # bmaj
            # Calculating the actual quantity
            targar[i][:,0] = (tar_bmaj_inter[i]+tar_bmaj_slope[i]*\
                              tar_bmaj_absc[i])*scale
            
            # Converting to pixel coordinates with dispersion
            targar[i][:,0] = targar[i][:,0]/targar[i][:,4]/np.sqrt(np.log(256))

            
            # bmin
            # Calculating the actual quantity
            targar[i][:,1] = (tar_bmin_inter[i]+tar_bmin_slope[i]*\
                              tar_bmin_absc[i])*scale
            
            # Converting to pixel coordinates with dispersion
            targar[i][:,1] = targar[i][:,1]/targar[i][:,4]/np.sqrt(np.log(256))
            
            # bpa
            # Calculating the actual quantity
            targar[i][:,2] = tar_bpa_inter[i]+tar_bpa_slope[i]*tar_bpa_absc[i]

            # Converting to rad
            targar[i][:,2] = np.pi*targar[i][:,2]/180.

            # Calculating sin
            targar[i][:,3] = np.sin(targar[i][:,2])

            # Calculating cos
            targar[i][:,4] = np.cos(targar[i][:,2])

        self._binfo_target = targar

        # Cascade down
        if self._gentrans_exe:
            self.gentrans(verb = verb)

        return

    def _initgentransvar(self, transnames = None, overwrite = None,
                         verb = None):
        """
        Check existence of variables, return True if a parameter is ill defined
        """
        output = False
        paras = locals().copy()
        paras.pop('self')
        for param in paras.keys():
            if type(paras[param]) != type(None):
                self.__dict__['_'+param] = copy.deepcopy(paras[param])
            else:
                if self.__dict__['_'+param] == None:
                    if verb or self._verb:
                        warnings.warn('Parameter {} is not defined.'.format(param))
                    output = True
        return output
            
    def gentrans(self, transnames = None, overwrite = None,
                  verb = None):
        """
        (De-)convolve input data cubes or images to target beam shapes

        Input:
        transnames (str or list of str): Output fits file names
        overwrite (bool)               : Overwrite output if already
                                         existent (True: yes)?

        Serially opens all cubes (images) listed in inputnames and
        generates cubes (De-)convolved to the resolution as listed in
        the target structure.

        """
        stop = self._initgentransvar(transnames = transnames,
                                     overwrite = overwrite, verb =
                                     verb)
        if stop:
            if verb or self._verb:
                warnings.warn('Parameters missing. Not generating output data'+ \
                              'sets.')
            return

        if type(self._target) == type(None):
            if verb or self._verb:
                warnings.warn('Attempting to generate target')
            self.gentarget(verb = verb)
        if type(self._target) == type(None):
            if verb or self._verb:
                warnings.warn('Failed to generate target properties.'+ \
                              'Not generating output.')
            return

        # Also check if transnames have the same type and number as inputnames
        if type(self._transnames) == type(''):
            transn = [self._transnames]
        else:
            transn = self._transnames

        if type(self._cubenames) == type(''):
            cuben = [self._cubenames]
        else:
            cuben = self._cubenames

        if len(transn) != len(inputn):
            if verb or self._verb:
                warnings.warn('Number of input names not matching number of'+ \
                              'output data sets. Not generating output data.')
            return

        # After this marginal verification, we continue
        for i in range(len(cuben)):
            incubus = fits.open(cuben[i])
            incubus[0].image.reshape((incubus[0].header['NAXIS1'],
                                      incubus[0].header['NAXIS2'],
                                      len(self._binfo_input[i][:,0])))
            
        return
        
def printcubeinfo(cubename):
    print('Header info ',cubename)
    hdu1= fits.open(cubename)
    for helement in hdu1[0].header.keys():
        if 'BMAJ' in helement:
            print(helement, hdu1[0].header[helement])
        if 'BMIN' in helement:
            print(helement, hdu1[0].header[helement])
        if 'BPA' in helement:
            print(helement, hdu1[0].header[helement])
        if 'CDELT2' in helement:
            print(helement, hdu1[0].header[helement])
        if 'FREQ' in helement:
            print(helement, hdu1[0].header[helement])
        if 'CELLSCAL' in helement:
            print(helement, hdu1[0].header[helement])
        if 'CTYPE3' in helement:
            print(helement, hdu1[0].header[helement])
    hdu1.close()

def printbeachconts(beach):
    print('Input cube 1')
    print('bmaj: ', beach.binfo_input[0][:,0])
    print('bmin: ', beach.binfo_input[0][:,1])
    print('bpa : ', beach.binfo_input[0][:,2])
    print('freq: ', beach.binfo_input[0][:,3])
    print('dx  : ', beach.binfo_input[0][:,4])
    print('Input cube 2')
    print('bmaj: ', beach.binfo_input[1][:,0])
    print('bmin: ', beach.binfo_input[1][:,1])
    print('bpa : ', beach.binfo_input[1][:,2])
    print('freq: ', beach.binfo_input[1][:,3])
    print('dx  : ', beach.binfo_input[1][:,4])
    print()
    print('Pixel cube 1')
    print('bmaj    : ', beach.binfo_pixel[0][:,0])
    print('bmin    : ', beach.binfo_pixel[0][:,1])
    print('cos bpa : ', beach.binfo_pixel[0][:,2])
    print('sin bpa : ', beach.binfo_pixel[0][:,3])
    print('freq    : ', beach.binfo_pixel[0][:,4])
    print('Pixel cube 2')
    print('bmaj    : ', beach.binfo_pixel[1][:,0])
    print('bmin    : ', beach.binfo_pixel[1][:,1])
    print('cos bpa : ', beach.binfo_pixel[1][:,2])
    print('sin bpa : ', beach.binfo_pixel[1][:,3])
    print('freq    : ', beach.binfo_pixel[1][:,4])
    print()
    print(beach._binfo_target)

if __name__ == '__main__':

    print('')
    print('#############')
    print('#############')
    print('#############')
    print('Test 1: Frequency axis')
    print('#############')
    printcubeinfo('p1_02.fits')
    printcubeinfo('p2_02.fits')
                  
    params = { 'cubes': ['p1_02.fits',
                         'p2_02.fits']
    }
    beach = Beach(cubenames = params['cubes'])
    print('#####')
    print('#####')
    print('#####')
    print('Got here')
    print('#####')
    print('#####')
    print('#####')
    beach.genbstats()
    beach.gentarget()
    printbeachconts(beach)

    print('')
    print('#############')
    print('#############')
    print('#############')
    print('Test 2: Velocity axis')
    print('#############')
    printcubeinfo('p1_04.fits')
    printcubeinfo('p2_04.fits')
    params = { 'cubes': ['p1_04.fits',
                     'p2_04.fits']
    }
    beach = Beach(cubenames = params['cubes'])
    beach.genbstats()
    beach.gentarget()
    printbeachconts(beach)


    print('')
    print('#############')
    print('#############')
    print('#############')
    print('Test 3: 2 axes')
    print('#############')
    printcubeinfo('p1_07.fits')
    printcubeinfo('p2_07.fits')
    params = { 'cubes': ['p1_07.fits',
                     'p2_07.fits']
    }
    beach = Beach(cubenames = params['cubes'])
    beach.genbstats()
    beach.gentarget()
    printbeachconts(beach)

    print('')
    print('#############')
    print('#############')
    print('#############')
    print('Test 4: CELLSCAL 1/F')
    print('#############')
    printcubeinfo('p1_05.fits')
    printcubeinfo('p2_05.fits')
    params = { 'cubes': ['p1_05.fits',
                     'p2_05.fits']
    }
    beach = Beach(cubenames = params['cubes'])
    beach.genbstats()
    beach.gentarget()
    printbeachconts(beach)
