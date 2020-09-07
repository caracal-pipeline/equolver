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
import pyfftw

#from bokeh.layouts import gridplot
#from bokeh.plotting import figure, output_file, show
#import bokeh.layouts as bokeh_layouts
import bokeh.plotting as bokeh_plotting
import bokeh.io as bokeh_io

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
                       scaling = 'frequency',
                       stype = ['stdev', 'commonbeam'],
                       sample = 'all', percents = 90,
                       tolerance = 0.1, nsamps = 200, epsilon = 0.0005,
                       tar_bmaj_inter = ['bmaj', 'frequency', 'maximum', 'total'],
                       tar_bmaj_slope = ['bmaj', 'frequency',  'stdev', 'total'],
                       tar_bmaj_absc = 0.0,
                       tar_bmin_inter = ['bmin', 'frequency',  'commonbeam', 'total'],
                       tar_bmin_slope = ['bmin', 'frequency',  'stdev', 'total'],
                       tar_bmin_absc = 0.0,
                       tar_bpa_inter =  ['bpa', 'frequency',  'commonbeam', 'total'],
                       tar_bpa_slope =  ['bpa', 'frequency',  'stdev', 'total'],
                       tar_bpa_absc =  0.0,
                       tar_scaling = 'frequency', genbstats_exe = True,
                       gentarget_exe = True,
                       gentrans_exe = True, tra_fitsnames = None,
                       tra_overwrite = False, 
                       tra_commonbeam = True, tra_indibeam = True,
                       threads = 1, verb = False):
        """Private instance variables:
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

        _binfo_input (list) : List of arrays of point spread func-
                                tion, in the order of headernames and
                                headers each array of size (chans, 6),
                                where first column bmaj in deg, second
                                column bmin in deg, third column bpa
                                in deg, fourth column frequency in Hz,
                                fifth column pixel size in deg, sixth
                                column reference frequency divided by
                                frequency, chans is the number of
                                channels.
        _binfo_pixel (list)  : _binfo_input converted into pixel
                               scaling using dispersion instead of HPBW
        _bstats (dict)       : Dictionary containing all statistics
        _binfo_target (list)  : Target beam properties

        """
        self._initvars()
        self._verb = verb
        self._threads = threads
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

        for para in ['parameter', 'scaling', 'stype', 'sample',
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
            
        for para in [ 'gentrans_exe', 'tra_fitsnames',
                      'tra_overwrite', 'tra_commonbeam',
                      'tra_indibeam']:
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

        self._tra_fitsnames = None
        self._tra_overwrite = None
        self._tra_commonbeam = None
        self._tra_indibeam = None

        self._threads = None
        self._verb = True
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
    def scaling(self):
        """
        Return a copy of parameter
        """
        return copy.deepcopy(self._scaling)

    @scaling.setter
    def scaling(self, value):
        """
        Set scaling
        """
        self._scaling = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @scaling.deleter
    def parameter(self, value):
        self._scaling = 'frequency'
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
        self._tar_bmaj_inter = ['bmaj', 'frequency', 'median', 'all']
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
        self._tar_bmaj_slope = ['bmaj', 'frequency', 'madstdev', 'all']
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
        self._tar_bmin_inter = ['bmaj', 'frequency', 'median', 'all']
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
        self._tar_bmin_slope = ['bmaj', 'frequency', 'madstdev', 'all']
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
        self._tar_bpa_inter = ['bpa', 'frequency', 'median', 'all']
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
        self._tar_bpa_slope = ['bpa', 'frequency', 'madstdev', 'all']
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
        self._tar_scaling = 'frequency'
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
    def gentrans_exe(self):
        """
        Return a copy of gentrans_exe
        """
        return self._gentrans_exe

    @gentrans_exe.setter
    def gentrans_exe(self, value):
        """
        Set gentrans_exe
        """
        self._gentrans_exe = value
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @gentrans_exe.deleter
    def gentrans_exe(self):
        self._gentrans_exe = False
        return

    @property
    def tra_fitsnames(self):
        """
        Return a copy of _tra_fitsnames
        """
        return self._tra_fitsnames

    @gentrans_exe.setter
    def gentrans_exe(self, value):
        """
        Set gentrans_exe
        """
        self._tra_fitsnames = copy.deepcopy(value)
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @gentrans_exe.deleter
    def tra_fitsnames(self):
        self._tra_fitsnames = None
        return

    @property
    def tra_overwrite(self):
        """
        Return a copy of verb
        """
        return self._tra_overwrite

    @tra_overwrite.setter
    def tra_overwrite(self, value):
        """
        Set tra_overwrite
        """
        self._tra_overwrite = value
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_overwrite.deleter
    def tra_overwrite(self):
        self._tra_overwrite = None
        return
    
    @property
    def tra_commonbeam(self):
        """
        Return a copy of verb
        """
        return self._tra_commonbeam

    @tra_commonbeam.setter
    def tra_commonbeam(self, value):
        """
        Set tra_commonbeam
        """
        self._tra_commonbeam = value
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_commonbeam.deleter
    def tra_commonbeam(self):
        self._tra_commonbeam = True
        return

    @property
    def tra_indibeam(self):
        """
        Return a copy of tra_indibeam
        """
        return self._tra_indibeam

    @tra_indibeam.setter
    def tra_indibeam(self, value):
        """
        Set tra_indibeam
        """
        self._tra_indibeam = value
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_indibeam.deleter
    def tra_indibeam(self):
        self._tra_indibeam = True
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

    @property
    def threads(self):
        """
        Return a copy of verb
        """
        return self._threads

    @threads.setter
    def threads(self, value):
        """
        Set threads
        """
        self._threads = threads
        return

    @threads.deleter
    def threads(self):
        self._threads = None
        return
    
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
        

            self._binfo_input.append(np.empty((naxis3,8,)))
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
                            restfreq = None, restfreq_replace = None,
                            normfreq = None, verb = False):
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
                    if self._verb or verb:
                        warnings.warn('Parameter {} is not defined.'.format(param))
                    output = True
        return output
    
    def genbinfo(self, bmaj = None, bmaj_replace = None,
                       bmin = None, bmin_replace = None,
                       bpa = None,  bpa_replace = None,
                       restfreq = HIFREQ, restfreq_replace = None,
                       normfreq = None, verb = True):
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
        if self._quantity is not None.
        """
        stop = self._initbinfo_inputvar(bmaj = bmaj,
                                        bmaj_replace = bmaj_replace,
                                        bmin = bmin, bmin_replace =
                                        bmin_replace, bpa = bpa,
                                        bpa_replace = bpa_replace,
                                        restfreq = restfreq,
                                        restfreq_replace =
                                        restfreq_replace,
                                        normfreq = normfreq, verb = verb)

        if stop:
            if self._verb:
                warnings.warn('Parameters missing. Not generating beam info.')
            return

        if type(self._headers) == type(None):
            if self._verb:
                warnings.warn('No headers loaded. Returning.')
            return
        
        if type(self._binfo_input) == type(None):
            self._initbinfo_input()
            
        if type(self._binfo_input) == type(None):
            if self._verb:
                warning('binfo_input not present. Use initheaders first.')
            return
        
        # Do not try to regenerate, if this is none, it has deliberately been
        # set to None, as __init__ has a default of 1E9
        if type(self._normfreq) == type(None):
            if verb or self._verb:
                warnings.warn('No normalization frequency read in, which '+ \
                              'disables further processing.')
            return

        normfreq = self._getdefault(self._normfreq)

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
            self._binfo_input[i][:,5] = np.ones(self._binfo_input[i].shape[0])
            if cellscal_use_constant:
                if not type(wcshand2) == type(None):
                    if not np.isnan(self._binfo_input[i][:,3].sum()):
                        self._binfo_input[i][:,5] = wcshand2.crval[2]/ \
                            self._binfo_input[i][:,3]

            # finlist: slowest index cube in 'headers' list, followed
            # by bmaj, bmin, and bpa
            self._binfo_input[i][:,0] = self._getchanval( \
                'BMAJ', self._headers[i], value = bmaj[i],
                 usevalue = bmaj_replace, usedefault = True, dscal = self._binfo_input[i][:,5])
            self._binfo_input[i][:,1] = self._getchanval( \
                'BMIN', self._headers[i], value = bmin[i],
                usevalue = bmin_replace, usedefault = True, dscal = self._binfo_input[i][:,5])
            self._binfo_input[i][:,2] = self._getchanval( \
                'BPA', self._headers[i], value = bpa[i],
                usevalue = bpa_replace, usedefault = True)

            if np.isnan(self._binfo_input[i][:,3].sum()):
                thafreq = normfreq[i]
            else:
                thafreq = self._binfo_input[i][:,3]
            self._binfo_input[i][:,6] = self._binfo_input[i][:,0]*\
                thafreq/normfreq[i]
            self._binfo_input[i][:,7] = self._binfo_input[i][:,1]*\
                thafreq/normfreq[i]
            
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

            # Scale to rad
            boutfarray[:,2] = np.pi*binfarray[:,2]/180.
            
            boutfarray[:,6] = binfarray[:,6]/binfarray[:,4]/np.sqrt(np.log(256))
            boutfarray[:,7] = binfarray[:,7]/binfarray[:,4]/np.sqrt(np.log(256))

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
            self.genbinfo(verb = verb)

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
            for key1 in ['constant', 'frequency']:
                self._bstats[key0][key1] = {}

                # stype
                for key2 in ['minimum', 'maximum', 'average', 'stdev', 'median',
                             'mad', 'madstdev', 'percentile', 'percents',
                             'commonbeam']:
                    self._bstats[key0][key1][key2] = {}
                    # for key4 in ['cube', 'chan', 'tota']:

                    # sample
                    self._bstats[key0][key1][key2]['total'] = np.nan
                    self._bstats[key0][key1][key2]['chan'] = \
                        np.empty((channels,))
                    self._bstats[key0][key1][key2]['chan'][:] = np.nan
                    self._bstats[key0][key1][key2]['cube'] = []
                    for i in range(ncubes):
                        self._bstats[key0][key1][key2]['cube'].append(np.nan)
                        
        return
    
    def _initbstatsvar(self, parameter = None, stype = None, scaling = None,
                       sample = None, percents = None, tolerance = None,
                       nsamps = None, epsilon = None, verb = False):
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
            
    def genbstats(self, parameter = None, scaling = None, stype = None,
                  sample = None, percents = None, tolerance = None,
                  nsamps = None, epsilon = None, verb = True):
        """Generate statistics and dump it into the bstats structure
        
        Input:
        parameter (str or list of str): Parameter name ('bmaj', 
                                        'bmin', 'bpa')
        scaling (str or list of str)    : Scaling type ('constant', 
                                        'frequency')
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
        the default), all fields are filled. If the scaling type is
        'constant', the given values for bmaj and bmin are evaluated,
        if 'frequency' is chosen, all beam sizes are scaled to the
        same norm-frequency nf, assuming that the beam sizes scale
        with 1/frequency. b(nf) = b(f)*f/nf .

        From top to bottom level:

        parameter: 
            major axis dispersions/hpbws ('bmaj')
            minor axis dispersions/hpbws ('bmin')
            beam position angles ('bpa')
        scaling:
            constant ('const')
            frequency ('frequency')
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

        stop = self._initbstatsvar(parameter = parameter, scaling =
                                   scaling, stype = stype, sample =
                                   sample, percents = percents,
                                   tolerance = tolerance, nsamps =
                                   nsamps, epsilon = epsilon, verb =
                                   verb)

        if stop:
            if self._verb:
                warnings.warn('Parameters missing. Not generating statistics.')
            return

        if self._bstats == None:
            if self._verb:
                warnings.warn('bstats not available, initializing bstats')
            self._initbstats(verb = self._verb)
        
        if self._bstats == None:
            if self._verb:
                warnings.warn('Failing to initialize bstats. Returning.')
            return
        
        para = copy.deepcopy(self._parameter)
        scal = copy.deepcopy(self._scaling)
        styp = copy.deepcopy(self._stype)
        samp = copy.deepcopy(self._sample)

        # Expand all parameters to the same format
        if type(para) == type(''):
            para = [para]
        if type(scal) == type(''):
            scal = [scal]
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
        if 'all' in scal:
            scal = ['constant', 'frequency']
        if 'all' in styp:
            styp = ['minimum', 'maximum', 'average', 'stdev', 'median',
                    'mad', 'madstdev', 'percentile', 'percents',
                    'commonbeam']
        if 'all' in samp:
            samp = ['cube', 'chan', 'total']
        
        # Make sure we do have the required information available
        if type(self._binfo_input) == type(None):
            if verb or self._verb:
                warnings.warn('Input information struct not '+ \
                              'available, regenerating.')
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
        for sca in scal:
            
            for sam in samp:

                # collect sample as list of np arrays
                collist = []

                #print(sam)

                if sam == 'cube':
                    for i in self._binfo_input:
                        if sca == 'constant':
                            collist.append(i[:,[0,1,2]])
                        else:
                            collist.append(i[:,[6,7,2]])
                            
                elif sam == 'chan':
                    for i in range(self._bstats['bmaj']['constant']['maximum'][sam].shape[0]):
                        collist.append(np.empty((0,3,), dtype=np.float64))
                        for j in self._binfo_input:
                            if j.shape[0] > i:
                                if sca == 'constant':
                                    collist[i] = np.append(collist[i],
                                                       j[i:i+1,[0,1,2]],
                                                       axis = 0)
                                else:
                                    collist[i] = np.append(collist[i],
                                                       j[i:i+1,[6,7,2]],
                                                       axis = 0)
                elif sam == 'total':
                    collist = [np.empty((0,3,))]
                    for i in self._binfo_input:
                        if sca == 'constant':
                            collist[0] = np.append(collist[0], i[:,[0,1,2]],
                                                   axis = 0)
                        else:
                            collist[0] = np.append(collist[0], i[:,[6,7,2]],
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
                        self._bstats[par][sca][sty][sam] = stats[:,parma[par]]
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
        
    def _getar(self, item, verb = False):
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
            if len(item) == 4:
                # The alternative is either numeric or quantity,
                # so this is sufficient for identifying the input
                if type(item[0]) == type(''):

                    # Ensure that we can actually do this
                    if type(self._bstats) == type(None) or \
                       np.isnan(self._bstats[item[0]][item[1]][item[2]][item[3]]):
                        if verb or self._verb:
                            warnings.warn('Parts of beam statistics struct '+ \
                                          'not available, regenerating.')
                            
                        # If the item is not present, we add it and
                        # generate the statistics
                        if type(self._parameter) == type(''):
                            if self._parameter != item[0]:
                                self._parameter = [self._parameter]
                        if type(self._scaling) == type(''):
                            if self._scaling != item[1]:
                                self._scaling = [self._scaling]
                        if type(self._stype) == type(''):
                            if self._stype != item[2]:
                                self._stype = [self._stype]
                        if type(self._sample) == type(''):
                            if self._sample != item[2]:
                                self._sample = [self._sample]
                        if not item[0] in self._parameter:
                            self._parameter += [item[0]]
                        if not item[1] in self._scaling:
                            self._scaling += [item[0]]
                        if not item[0] in self._stype:
                            self._stype += [item[0]]
                        if not item[0] in self._sample:
                            self._sample += [item[0]]
                        self.genbstats(parameter = item[0], scaling = item[1], stype =
                                       item[2], sample = item[3], verb = verb)
                    if type(self._bstats) == type(None) or \
                       np.isnan(self._bstats[item[0]][item[1]][item[2]][item[3]]):
                        if verb or self._verb:
                            warnings.warn('Failed to genereate beam statistics'+ \
                                          'struct. Returning empty-handed.')
                        return               

                    # This should have cascaded up, such that the information
                    # is at our hands
                    output = []
                    for i in range(len(self._headers)):
                        output.append(np.zeros(self._binfo_input[i].shape[0]))

                        if item[3] == 'cube':
                            output[i] = output[i]+\
                                self._bstats[item[0]][item[1]][item[2]][item[3]][i]
                        elif item[3] == 'channel' or item[3] == 'total':
                            output[i] = output[i]+\
                                self._bstats[item[0]][item[1]][item[2]][item[3]]
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
                               the target array, either 'frequency' or
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
                print('er')
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
                warnings.warn('bstats not available, regenerating.')
            self.genbstats(verb = self._verb)
        if type(self._bstats) == type(None):
            if verb or self._verb:
                warnings.warn('Failed to generate bstats. Not generating target.')
            return
        tar_bmaj_inter = self._getar(self._tar_bmaj_inter, verb = self._verb)
        tar_bmaj_slope = self._getar(self._tar_bmaj_slope, verb = self._verb)
        tar_bmaj_absc  = self._getar(self._tar_bmaj_absc, verb = self._verb)
        tar_bmin_inter = self._getar(self._tar_bmin_inter, verb = self._verb)
        tar_bmin_slope = self._getar(self._tar_bmin_slope, verb = self._verb)
        tar_bmin_absc  = self._getar(self._tar_bmin_absc, verb = self._verb)
        tar_bpa_inter  = self._getar(self._tar_bpa_inter, verb = self._verb)
        tar_bpa_slope  = self._getar(self._tar_bpa_slope, verb = self._verb)
        tar_bpa_absc   = self._getar(self._tar_bpa_absc, verb = self._verb)

        # For the preliminary output we make a copy of input
        targar = copy.deepcopy(self._binfo_pixel)

        for i in range(len(targar)):

            # Notice that this is the inverse operation to translating
            # the input beams to frequency scaling
            if np.isnan(self._binfo_input[i][:,3].sum()):
                thafreq = normfreq[i]
            else:
                thafreq = targar[i][:,3]

            if self._tar_scaling == 'frequency':
                scale = normfreq[i]/thafreq
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

        self._binfo_target = targar

        # Cascade down
        if self._gentrans_exe:
            self.gentrans(verb = verb)

        return

    def _initgentransvar(self, tra_fitsnames = None, tra_overwrite =
                         None, tra_commonbeam = None, tra_indibeam =
                         None, threads = None, verb = None):
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
            
    def gentrans(self, tra_fitsnames = None, tra_overwrite = None,
                       tra_commonbeam = None, tra_indibeam = None,
                       threads = None, verb = True):
        """
        (De-)convolve input data cubes or images to target beam shapes

        Input:
        tra_fitsnames (str or list of str): Output fits file names
        tra_overwrite (bool)               : Overwrite output if already
                                         existent (True: yes)?

        Serially opens all cubes (images) listed in inputnames and
        generates cubes (De-)convolved to the resolution as listed in
        the target structure.

        """
        stop = self._initgentransvar(tra_fitsnames = tra_fitsnames,
                                     tra_overwrite = tra_overwrite,
                                     tra_commonbeam = tra_commonbeam,
                                     tra_indibeam = tra_indibeam, verb
                                     = verb, threads = threads)
        if stop:
            if verb or self._verb:
                warnings.warn('Parameters missing. Not generating output data'+ \
                              'sets.')
            return

        if type(self._binfo_target) == type(None):
            if verb or self._verb:
                warnings.warn('Target information not available, regenerating.')
            self.gentarget(verb = self._verb)
        if type(self._binfo_target) == type(None):
            if verb or self._verb:
                warnings.warn('Failed to generate target properties.'+ \
                              'Not generating output.')
            return

        # Also check if tra_fitsnames have the same type and number as inputnames
        if type(self._tra_fitsnames) == type(''):
            transn = [self._tra_fitsnames]
        else:
            transn = self._tra_fitsnames

        if type(self._cubenames) == type(''):
            cuben = [self._cubenames]
        else:
            cuben = self._cubenames

        if len(transn) != len(cuben):
            if verb or self._verb:
                warnings.warn('Number of input names not matching number of'+ \
                              'output data sets. Not generating output data.')
            return

        # After this marginal verification, we continue

        # Open, reconvolve, copy
        for i in range(len(cuben)):
            incubus = fits.open(cuben[i])
            incubus_image = incubus[0].data.astype('float'+'{:d}'.format(incubus[0].data.itemsize*8))
            orishape = incubus_image.shape
            #incubus_image = incubus_image.reshape((incubus[0].header['NAXIS1'],
            #                             incubus[0].header['NAXIS2'],
            #                             len(self._binfo_input[i][:,0])))
            if incubus_image.ndim > 3:
                incubus_image = np.squeeze(incubus_image)
            if incubus_image.ndim == 2:
                incubus_image = incubus_image.reshape((1, incubus_image.shape[0], incubus_image.shape[1]))
            # Make a copy
            outcubus_image = incubus_image.copy()*0.+np.nan
            
            for plane in range(incubus_image.shape[0]):
                print('Processing {:s} plane {:d}'.format(cuben[i],plane))
                originbeam = self._binfo_pixel[i][plane,:]
                targetbeam = self._binfo_target[i][plane,:]
                originplane = incubus_image[plane,:,:]
                targetplane = outcubus_image[plane,:,:]
                self._reconvolve(originplane, targetplane, originbeam,
                                 targetbeam, threads = self._threads)

            # Replace incubus with outcubus and write out
            incubus[0].data[:] = outcubus_image.astype(incubus[0].data.dtype)

            if self._tra_commonbeam:
                if self._tar_scaling == 'frequency':

                    # For each cube calculate the average beam
                    # normalised to the reference frequency dscal is
                    # the scale per cube w.r.t. the reference pixel,
                    # reference frequency divided by frequency for
                    # each channel
                    divisor = self._binfo_input[i][:,5]
                    beamscal = '1/F'
                else:
                    divisor = 1.
                    beamscal = 'CONSTANT'

                # Divide by reference frequency divided by
                # frequency, then multiply with pixel size, then
                # take average, if user intelligent, then the average
                # is identical to the scaled value for every channel
                bmajav = np.sqrt(np.log(256))*np.average((self._binfo_target[i][:,0]/divisor)*
                                    self._binfo_target[i][:,4])
                bminav = np.sqrt(np.log(256))*np.average((self._binfo_target[i][:,1]/divisor)*
                                    self._binfo_target[i][:,4])
                bpav = 180.*np.average(self._binfo_target[i][:,2])/np.pi
                
                incubus[0].header['BMAJ'] = bmajav
                incubus[0].header['BMIN'] = bminav
                incubus[0].header['BPA'] = bpav
                incubus[0].header['BEAMSCAL'] = beamscal
            
            if self._tra_indibeam:
                
                # Here just copy all target beam properties into the header
                for j in range(len(self._binfo_input[i][:,0])):
                    incubus[0].header['BMAJ{:d}'.format(j+1)] = np.sqrt(np.log(256))*self._binfo_target[i][j,0]*\
                                    self._binfo_target[i][j,4]
                    incubus[0].header['BMIN{:d}'.format(j+1)] = np.sqrt(np.log(256))*self._binfo_target[i][j,1]*\
                                    self._binfo_target[i][j,4]
                    incubus[0].header['BPA{:d}'.format(j+1)] = 180.*self._binfo_target[i][j,2]/np.pi

            # Finally write and close cube
            incubus.writeto(self._tra_fitsnames[i], overwrite = self._tra_overwrite)
            incubus.close()
        return

    def _gaussian_2dp(self, naxis1 = 100, naxis2 = 100, cdelt1 = 1.,
                      cdelt2 = 1., amplitude_maj_a = 1.,
                      dispersion_maj_a = np.inf, signum_maj_a = -1.,
                      amplitude_min_a = 1., dispersion_min_a = np.inf,
                      signum_min_a = -1., pa_a = 0., amplitude_maj_b =
                      1., dispersion_maj_b = np.inf, amplitude_min_b =
                      1., dispersion_min_b = np.inf, signum_maj_b =
                      1., signum_min_b = 1., pa_b = 0., dtype =
                      'float32', centering = 'origin', forreal = True):
        """
        Returns the the product of two Quasi-Gaussians as ndarray
        (positve sign in exponent allowed)

        Input:
        naxis1 (int)            : Number of pixels axis 1 (FITS 
                                  convention)
        naxis2 (int)            : Number of pixels axis 2 (FITS
                                  convention)
        cdelt1 (float)          : Pixel size axis 1 (FITS convention)
        cdelt2 (float)          : Pixel size axis 2 (FITS convention)
        amplitude_maj_a (float) : Amplitude Gaussian a, major axis, 
                                  should be > 0, can be np.inf
        dispersion_maj_a (float): Dispersion Gaussian a, major axis
        signum_maj_a (int)      : Signum in the exponent for major axis
                                  component
        amplitude_min_a (float) : Amplitude Gaussian a, minor axis,
                                  should be > 0, can be np.inf
        dispersion_min_a (float): Dispersion Gaussian a, minor axis
        signum_min_a (int)      : Signum in the exponent for minor axis 
                                  component
        pa_a (float)            : Position angle of half major axis of
                                  Gaussian a, measured from positive 
                                  direction of axis 2 through origin to
                                  half major axis of Gaussian in rad
        amplitude_maj_b (float) : Amplitude Gaussian b, major axis,
                                  should be > 0, can be np.inf
        dispersion_maj_b (float): Dispersion Gaussian b, major axis
        signum_maj_b (int)      : Signum in the exponent for major axis
                                  component
        amplitude_min_b (float) : Amplitude Gaussian b, minor axis,
                                  should be > 0, can be np.inf
        dispersion_min_b (float): Dispersion Gaussian b, minor axis
        signum_min_b (int)      : Signum in the exponent for minor axis 
                                  component
        pa_b (float)            : Position angle of half major axis of
                                  Gaussian b, measured from positive 
                                  direction of axis 2 through origin to
                                  half major axis of Gaussian in rad
        dtype (str)             : numpy dtype 
        centering = 'origin'    : centre on 'origin' (pixel 0,0) and
                                  assume that for pixel > naxisi//2
                                  pixel = pixel-naxisi or on 'centre'
                                  to place the Gaussian at naxisi//2,
                                  alternatively provide center as a 
                                  pair of float [x1,x2]
        forreal                 : Only relevant for centering = 
                                  'origin'. Assume the map to be the 
                                  result of a real Fourier 
                                  transformation and hence do not make
                                  the 'origin' assumption for axis 1

        Calculates the product of four planar quasi-Gaussians

            P = amplitude_maj_a*amplitude_maj_b*
                signum_maj_a/2*(x0a/dispersion_maj_a)+
                signum_min_a/2*(x1a/dispersion_min_a)+
                amplitude_min_a*amplitude_min_b*
                signum_maj_b/2*(x0b/dispersion_maj_b)+
                signum_min_b/2*(x1b/dispersion_min_b)
                
        where
            x1a =  cos(pa_a)*x1+sin(pa_a)*x0
            x0a = -sin(pa_a)*x1+cos(pa_a)*x0
            x1b =  cos(pa_b)*x1+sin(pa_b)*x0
            x0b = -sin(pa_b)*x1+cos(pa_b)*x0

        (Notice that x0 is measured along the conventional y-axis or
        naxis2 and x1 along the x-axis of naxis1). In case of
        dispersion_xxx_y == np.inf, a constant is the result. If
        centering == "origin", the result is centered on the origin,
        which is the correct approach if the function is supposed to
        be (the Fourier transform of) a convolution kernel.

        """

        # Create map 
        indices = np.indices((naxis2,naxis1), dtype = dtype)

        if centering == 'origin':
            
            # Zero at origin and negative half way through
            indices[0][indices[0] > naxis2/2] = indices[0][indices[0] > naxis2/2]-naxis2
            if not forreal:

                # If for a real transformation only one axis
                indices[1][indices[1] > naxis1/2] = indices[1][indices[1] > naxis1/2]-naxis1
                
        elif type(centering) == type([]):
            indices[0] = indices[0] - centering[1]
            indices[1] = indices[1] - centering[0]          
        else:
            # Centered on the centre of the map
            indices[0] = indices[0] - naxis2//2.
            indices[1] = indices[1] - naxis1//2.

        # Scaling
        indices[0][:] = cdelt2*indices[0][:]
        indices[1][:] = cdelt1*indices[1][:]
        
        xmina =  np.cos(pa_a)*indices[1]+np.sin(pa_a)*indices[0]
        xmaja = -np.sin(pa_a)*indices[1]+np.cos(pa_a)*indices[0]
        xminb =  np.cos(pa_b)*indices[1]+np.sin(pa_b)*indices[0]
        xmajb = -np.sin(pa_b)*indices[1]+np.cos(pa_b)*indices[0]

        np.seterr(divide = 'ignore')
        exponent1a = signum_maj_a*np.power(np.divide(xmaja,dispersion_maj_a),2)
        exponent0a = signum_min_a*np.power(np.divide(xmina,dispersion_min_a),2)
        exponent1b = signum_maj_b*np.power(np.divide(xmajb,dispersion_maj_b),2)
        exponent0b = signum_min_b*np.power(np.divide(xminb,dispersion_min_b),2)

        np.seterr(divide = None)

        amplitude = amplitude_maj_a*amplitude_min_a*amplitude_maj_b*amplitude_min_b
        exponent=exponent1a+exponent0a+exponent1b+exponent0b
        
        return amplitude*np.exp(0.5*(exponent1a+exponent0a+exponent1b+exponent0b))

    def _igaussian_2dp(self, naxis1 = 100, naxis2 = 100, cdelt1 = 1.,
                      cdelt2 = 1., amplitude_maj_a = 1.,
                      dispersion_maj_a = 0., signum_maj_a = -1.,
                      amplitude_min_a = 1., dispersion_min_a = 0.,
                      signum_min_a = -1., pa_a = 0., amplitude_maj_b =
                      1., dispersion_maj_b = 0., amplitude_min_b =
                      1., dispersion_min_b = 0., signum_maj_b =
                      1., signum_min_b = 1., pa_b = 0., dtype =
                      'float32', centering = 'origin', forreal = True):
        """Returns the the coefficients of the inverst FT of a product of two Quasi-Gaussians as ndarray
        (positve sign in exponent allowed)

        Same as _gaussian_2dp, only assuming that the Fourier
        transform instead of the original Gaussians are
        constructed. This means that cdelt and the dispersions are
        adjusted. The Array returned can be used in Fourier space

        """

        cdelt1 = 1./(cdelt1*naxis1)
        cdelt2 = 1./(cdelt2*naxis2)
        
        if dispersion_maj_a != 0.:
            amplitude_maj_a = amplitude_maj_a*dispersion_maj_a*np.sqrt(2.*np.pi)
        if dispersion_min_a != 0.:
            amplitude_min_a = amplitude_min_a*dispersion_min_a*np.sqrt(2.*np.pi)
        if dispersion_maj_b != 0.:
            amplitude_maj_b = amplitude_maj_b*dispersion_maj_b*np.sqrt(2.*np.pi)
        if dispersion_min_b != 0.:
            amplitude_min_b = amplitude_min_b*dispersion_min_b*np.sqrt(2.*np.pi)

        amplitude_maj_a = np.power(amplitude_maj_a,-signum_maj_a)
        amplitude_min_a = np.power(amplitude_min_a,-signum_min_a)

        amplitude_maj_b = np.power(amplitude_maj_b,-signum_maj_b)
        amplitude_min_b = np.power(amplitude_min_b,-signum_min_b)
               
        # We expressively allow for zero dispersion
        np.seterr(divide = 'ignore')
        if dispersion_maj_a == 0.:
            dispersion_maj_a = np.PZERO
        if dispersion_maj_b == 0.:
            dispersion_maj_b = np.PZERO
        if dispersion_min_a == 0.:
            dispersion_min_a = np.PZERO
        if dispersion_min_b == 0.:
            dispersion_min_b = np.PZERO

        dispersion_maj_a = np.divide(1.,2.*np.pi*dispersion_maj_a)
        dispersion_min_a = np.divide(1.,2.*np.pi*dispersion_min_a)
        dispersion_maj_b = np.divide(1.,2.*np.pi*dispersion_maj_b)
        dispersion_min_b = np.divide(1.,2.*np.pi*dispersion_min_b)
        np.seterr(divide = None)

        if forreal:
            naxis1 = naxis1//2+1

        return self._gaussian_2dp(naxis1 = naxis1, naxis2 = naxis2,
                                  cdelt1 = cdelt1, cdelt2 = cdelt2,
                                  amplitude_maj_a = amplitude_maj_a,
                                  dispersion_maj_a = dispersion_maj_a,
                                  signum_maj_a = signum_maj_a,
                                  amplitude_min_a = amplitude_min_a,
                                  dispersion_min_a = dispersion_min_a,
                                  signum_min_a = signum_min_a,
                                  pa_a = pa_a,
                                  amplitude_maj_b = amplitude_maj_b,
                                  dispersion_maj_b = dispersion_maj_b,
                                  amplitude_min_b = amplitude_min_b,
                                  dispersion_min_b = dispersion_min_b,
                                  signum_maj_b = signum_maj_b, signum_min_b
                                  = signum_min_b, pa_b = pa_b, dtype =
                                  dtype, centering = centering, forreal =
                                  forreal)

    def convoltests(self, point_source = 'point_source.fits', gaussian_at_centre = 'gaussian_at_centre.fits', gaussian_at_origin = 'gaussian_at_origin.fits', real_fft_conv = 'real_fft_conv.fits', real_fft_conv_calc = 'real_fft_conv_calc.fits', reconvolve_input_image = 'reconvolve_input_image.fits', reconvolve_output_image = 'reconvolve_output_image.fits'):

        threads = 1
        
        print()
        print('###############################')
        print('###############################')
        print('###############################')
        print(' Exercises in FFT convolutions')
        print('###############################')
        print()
        print('###############################')
        print('Make a map with Point source at centre')
        print('###############################')
        print()

        newar = np.zeros((1025,513), dtype = '>f4')
        newar[newar.shape[0]//2,newar.shape[1]//2] = 1.

        hdu = fits.HDUList([fits.PrimaryHDU(newar)])
        target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))

        hdu[0].data[:] = target.astype(newar.dtype)
        print('Result at centre', hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2])

        hdu.writeto(point_source, overwrite = True)
        print('Image to be found at', point_source)
        hdu.close()

        print()
        print('#########################################')
        print('Make a map with Gaussian at centre')
        print('#########################################')
        print()

        newar = np.zeros((1025,513), dtype = '>f4')
        hdu = fits.HDUList([fits.PrimaryHDU(newar)])
        target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))

        # Generate Gaussian
        sign_maj_a = -1.
        amp_maj_a = 1.
        HPBW_maj_a = 8.
        sign_min_a = -1.
        amp_min_a = 1.
        HPBW_min_a = 8.
        pang_a = 30.
        
        dis_maj_a = HPBW_maj_a/np.sqrt(np.log(256.))
        dis_min_a = HPBW_min_a/np.sqrt(np.log(256.))
        pa_a = np.pi*30./180.

        hdu[0].data[:] = self._gaussian_2dp( naxis1 =
                                             target.shape[1], naxis2 =
                                             target.shape[0], cdelt1 =
                                             1., cdelt2 = 1.,
                                             amplitude_maj_a =
                                             amp_maj_a,
                                             dispersion_maj_a =
                                             dis_maj_a, signum_maj_a =
                                             sign_maj_a,
                                             amplitude_min_a =
                                             amp_min_a,
                                             dispersion_min_a =
                                             dis_min_a, signum_min_a =
                                             sign_min_a, pa_a = pa_a,
                                             dtype = target.dtype,
                                             centering = 'bla',
                                             forreal = False).astype(hdu[0].data.dtype)
        
        print('result at centre', hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2])

        xmin = np.cos(pa_a)*HPBW_min_a/2
        ymin = np.sin(pa_a)*HPBW_min_a/2
        xmaj = -np.sin(pa_a)*HPBW_maj_a/2
        ymaj = np.cos(pa_a)*HPBW_maj_a/2

        print('result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)])
        print('result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)])
        
        hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2] += 2
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)] += 2.
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)] += 2.
        hdu.writeto(gaussian_at_centre, overwrite = True)
        print('Image to be found at', gaussian_at_centre)
        hdu.close()
        
        print()
        print('#########################################')
        print('Generate Gaussian at origin')
        print('#########################################')
        print()
        
        newar = np.zeros((1025,513), dtype = '>f4')
        hdu = fits.HDUList([fits.PrimaryHDU(newar)])
        target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))

        # Generate Gaussian
        sign_maj_a = -1.
        amp_maj_a = 1.
        HPBW_maj_a = 8.
        sign_min_a = -1.
        amp_min_a = 1.
        HPBW_min_a = 8.
        pang_a = 30.
        
        dis_maj_a = HPBW_maj_a/np.sqrt(np.log(256.))
        dis_min_a = HPBW_min_a/np.sqrt(np.log(256.))
        pa_a = np.pi*30./180.

        hdu[0].data[:] = self._gaussian_2dp( naxis1 =
                                             target.shape[1], naxis2 =
                                             target.shape[0], cdelt1 =
                                             1., cdelt2 = 1.,
                                             amplitude_maj_a =
                                             amp_maj_a,
                                             dispersion_maj_a =
                                             dis_maj_a, signum_maj_a =
                                             sign_maj_a,
                                             amplitude_min_a =
                                             amp_min_a,
                                             dispersion_min_a =
                                             dis_min_a, signum_min_a =
                                             sign_min_a, pa_a = pa_a,
                                             dtype = target.dtype,
                                             centering = 'origin',
                                             forreal = False).astype(hdu[0].data.dtype)
        
        print('result at centre', hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2])

        xmin = np.cos(pa_a)*HPBW_min_a/2
        ymin = np.sin(pa_a)*HPBW_min_a/2
        xmaj = -np.sin(pa_a)*HPBW_maj_a/2
        ymaj = np.cos(pa_a)*HPBW_maj_a/2

        print('result at half power', hdu[0].data[int(ymin),int(xmin)])
        print('result at half power', hdu[0].data[int(ymaj),int(xmaj)])
        
        hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2] += 2.
        hdu[0].data[int(ymin),int(xmin)] += 2.
        hdu[0].data[int(ymaj),int(xmaj)] += 2.
        hdu.writeto(gaussian_at_origin, overwrite = True)
        print('Image to be found at', gaussian_at_origin)
        hdu.close()
        print('')
        
        print('#########################################')
        print('Do a real FFT convolution')
        print('#########################################')
        print()

        newar = np.zeros((1025,513), dtype = '>f4')
        newar[newar.shape[0]//2,newar.shape[1]//2] = 1.
        
        hdu = fits.HDUList([fits.PrimaryHDU(newar)])
        target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))

        # Generate Gaussian
        sign_maj_a = -1.
        amp_maj_a = 1.
        HPBW_maj_a = 8.
        sign_min_a = -1.
        amp_min_a = 1.
        HPBW_min_a = 8.
        pang_a = 30.
        
        dis_maj_a = HPBW_maj_a/np.sqrt(np.log(256.))
        dis_min_a = HPBW_min_a/np.sqrt(np.log(256.))
        pa_a = np.pi*30./180.

        kernel = self._gaussian_2dp(naxis1 = target.shape[1], naxis2 =
                                    target.shape[0], cdelt1 =
                                    1., cdelt2 = 1.,
                                    amplitude_maj_a =
                                    amp_maj_a,
                                    dispersion_maj_a =
                                    dis_maj_a, signum_maj_a =
                                    sign_maj_a,
                                    amplitude_min_a =
                                    amp_min_a,
                                    dispersion_min_a =
                                    dis_min_a, signum_min_a =
                                    sign_min_a, pa_a = pa_a,
                                    dtype = target.dtype,
                                    centering = 'origin',
                                    forreal = False)
        
        fft = pyfftw.builders.rfft2(kernel, planner_effort=None, threads=threads, auto_align_input=True, auto_contiguous=True, avoid_copy=False, norm=None)
        ikernel = fft()

        itarget = ikernel.copy()
        fft.update_arrays(target.astype(fft.input_dtype), itarget)
        fft()
        iconvolved = ikernel*itarget

        ifft = pyfftw.builders.irfft2(iconvolved, s=target.shape, planner_effort=None, threads=threads, auto_align_input=True, auto_contiguous=True, avoid_copy=False, norm=None)
        convolved = ifft()
        
        hdu[0].data[:] = convolved.astype(hdu[0].data.dtype)
        
        xmin = np.cos(pa_a)*HPBW_min_a/2
        ymin = np.sin(pa_a)*HPBW_min_a/2
        xmaj = -np.sin(pa_a)*HPBW_maj_a/2
        ymaj = np.cos(pa_a)*HPBW_maj_a/2

        print('result at centre', hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2])
        print('result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)])
        print('result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)])
        hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2] += 2
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)] += 2.
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)] += 2.
        hdu.writeto(real_fft_conv, overwrite = True)
        print('Image to be found at', real_fft_conv)
        hdu.close()
        print()
        
        print('#########################################')
        print('Do a real FFT convolution with a Gaussian calculated in the Fourier domain')
        print('#########################################')
        print()

        newar = np.zeros((1025,513), dtype = '>f4')
        newar[newar.shape[0]//2,newar.shape[1]//2] = 1.

        hdu = fits.HDUList([fits.PrimaryHDU(newar)])
        target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))

        # Generate Gaussian
        sign_maj_a = -1.
        amp_maj_a = 1.
        HPBW_maj_a = 8.
        sign_min_a = -1.
        amp_min_a = 1.
        HPBW_min_a = 8.
        pang_a = 30.
        
        dis_maj_a = HPBW_maj_a/np.sqrt(np.log(256.))
        dis_min_a = HPBW_min_a/np.sqrt(np.log(256.))
        pa_a = np.pi*30./180.

        ikernel = self._igaussian_2dp(naxis1 = target.shape[1], naxis2
                                      = target.shape[0], cdelt1 = 1.,
                                      cdelt2 = 1., amplitude_maj_a = amp_maj_a,
                                      dispersion_maj_a = dis_maj_a,
                                      signum_maj_a = sign_maj_a,
                                      amplitude_min_a = amp_min_a,
                                      dispersion_min_a = dis_min_a,
                                      signum_min_a = sign_min_a, pa_a = pa_a,
                                      dtype =
                                      target.dtype, centering =
                                      'origin', forreal = True)

        fft = pyfftw.builders.rfft2(target.copy(), planner_effort=None, threads=threads, auto_align_input=True, auto_contiguous=True, avoid_copy=False, norm=None)
        itarget = fft()

        iconvolved = ikernel*itarget

        ifft = pyfftw.builders.irfft2(iconvolved, s=target.shape, planner_effort=None, threads=threads, auto_align_input=True, auto_contiguous=True, avoid_copy=False, norm=None)
        convolved = ifft()
        
        hdu[0].data[:] = convolved.astype(hdu[0].data.dtype)
        
        print('result at centre', hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2])

        xmin = np.cos(pa_a)*HPBW_min_a/2
        ymin = np.sin(pa_a)*HPBW_min_a/2
        xmaj = -np.sin(pa_a)*HPBW_maj_a/2
        ymaj = np.cos(pa_a)*HPBW_maj_a/2
        
        print('result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)])
        print('result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)])
        hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2] += 2
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)] += 2.
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)] += 2.
        hdu.writeto(real_fft_conv_calc, overwrite = True)
        print('Image to be found at', real_fft_conv_calc)
        hdu.close()
        
        print()
        print('#########################################')
        print('Finally do the magic to turn an existing Gaussian into another')
        print('#########################################')
        print()

        newar = np.zeros((1025,513), dtype = '>f4')
        hdu = fits.HDUList([fits.PrimaryHDU(newar)])
        target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))

        # Generate Gaussian
        sign_maj_a = -1.
        amp_maj_a = 1.
        HPBW_maj_a = 8.
        sign_min_a = -1.
        amp_min_a = 1.
        HPBW_min_a = 8.
        pang_a = 30.
        
        dis_maj_a = HPBW_maj_a/np.sqrt(np.log(256.))
        dis_min_a = HPBW_min_a/np.sqrt(np.log(256.))
        pa_a = np.pi*30./180.

        hdu[0].data[:] = self._gaussian_2dp( naxis1 =
                                             target.shape[1], naxis2 =
                                             target.shape[0], cdelt1 =
                                             1., cdelt2 = 1.,
                                             amplitude_maj_a =
                                             amp_maj_a,
                                             dispersion_maj_a =
                                             dis_maj_a, signum_maj_a =
                                             sign_maj_a,
                                             amplitude_min_a =
                                             amp_min_a,
                                             dispersion_min_a =
                                             dis_min_a, signum_min_a =
                                             sign_min_a, pa_a = pa_a,
                                             dtype = target.dtype,
                                             centering = 'bla',
                                             forreal = False).astype(hdu[0].data.dtype)
        
        print('Original result at centre', hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2])

        xmin = np.cos(pa_a)*HPBW_min_a/2
        ymin = np.sin(pa_a)*HPBW_min_a/2
        xmaj = -np.sin(pa_a)*HPBW_maj_a/2
        ymaj = np.cos(pa_a)*HPBW_maj_a/2

        print('Original result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)])
        print('Original result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)])
        
        hdu.writeto(reconvolve_input_image, overwrite = True)
        print('Image pure Gaussian to be found at', reconvolve_input_image)
        hdu.close()


        hdu = fits.open(reconvolve_input_image)
        target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))

        # This is repeating the stuff above with an important
        # difference: signum of a is inverted
        sign_maj_a = -1.
        amp_maj_a = 1.
        HPBW_maj_a = 8.
        sign_min_a = -1.
        amp_min_a = 1.
        HPBW_min_a = 8.
        pang_a = 30.

        # Happens here
        sign_maj_a = -sign_maj_a
        sign_min_a = -sign_min_a

        dis_maj_a = HPBW_maj_a/np.sqrt(np.log(256.))
        dis_min_a = HPBW_min_a/np.sqrt(np.log(256.))
        pa_a = np.pi*pang_a/180.

        # This is now the Gaussian that we want
        sign_maj_b = -1.
        a_maj_b = 1.
        HPBW_maj_b = 9.
        sign_min_b = -1.
        a_min_b = 1.
        HPBW_min_b = 8.
        pang_b = 0.
        
        dis_maj_b = HPBW_maj_b/np.sqrt(np.log(256.))
        dis_min_b = HPBW_min_b/np.sqrt(np.log(256.))
        pa_b = np.pi*pang_b/180.

        # The following is an old strategy, (de-)convolving the
        # original image with the minimally required beam. Has been
        # moderately successful as de-convolution still allows for a
        # few pixels only for very small kernels. The problem with
        # this strategy so far is the unknown normalization.
        
        #minikern = np.amin([dis_min_a, dis_maj_a, dis_min_b, dis_maj_b])
        #disn_min_a = np.power(dis_min_a,2.)-np.power(minikern,2.)
        #disn_maj_a = np.power(dis_maj_a,2.)-np.power(minikern,2.)
        #disn_min_b = np.power(minikern,2.)-np.power(dis_min_b,2.)
        #disn_maj_b = np.power(minikern,2.)-np.power(dis_maj_b,2.)
        #
        #sign_min_a = 1.
        #sign_maj_a = 1.
        #sign_min_b = -1.
        #sign_maj_b = -1.
        #
        #if disn_min_a != 0.:
        #    sign_min_a = np.sign(disn_min_a)
        #if disn_maj_a != 0.:
        #    sign_maj_a = np.sign(disn_maj_a)
        #if disn_min_b != 0.:
        #    sign_min_b = np.sign(disn_min_b)
        #if disn_maj_b != 0.:
        #    sign_maj_b = np.sign(disn_maj_b)

        #disn_min_a = np.sqrt(np.abs(disn_min_a))
        #disn_maj_a = np.sqrt(np.abs(disn_maj_a))
        #disn_min_b = np.sqrt(np.abs(disn_min_b))
        #disn_maj_b = np.sqrt(np.abs(disn_maj_b))
        
        #bcorrmaj = np.power(np.sqrt(2*np.pi*disn_maj_b*disn_maj_b*disn_maj_a*disn_maj_a/(disn_maj_b*disn_maj_b+disn_maj_a*disn_maj_a)), sign_maj_b)
        #bcorrmin = np.power(np.sqrt(2*np.pi*disn_min_b*disn_min_b*disn_min_a*disn_min_a/(disn_min_b*disn_min_b+disn_min_a*disn_maj_a)), sign_min_b)
        #bcorrmin = np.power(np.sqrt(2*np.pi*disn_min_a*disn_min_a*minikern*minikern/(disn_min_b*disn_min_b+minikern*minikern)), sign_min_b)
        #acorrmaj = np.power(np.sqrt(2*np.pi*disn_maj_a*disn_maj_a*minikern*minikern/(disn_maj_a*disn_maj_a+minikern*minikern)), -1)
        #acorrmin = np.power(np.sqrt(2*np.pi*disn_min_a*disn_min_a*minikern*minikern/(disn_min_a*disn_min_a+minikern*minikern)), sign_min_a)

        #ikernel = self._igaussian_2dp(naxis1 = target.shape[1], naxis2
        #                             = target.shape[0], cdelt1 = 1.,
        #                             cdelt2 = 1., amplitude_maj_a =
        #                             1., dispersion_maj_a =
        #                             disn_maj_a, signum_maj_a =
        #                             sign_maj_a, amplitude_min_a = 1.,
        #                             dispersion_min_a = disn_min_a,
        #                             signum_min_a = sign_min_a, pa_a =
        #                             pa_a, amplitude_maj_b = 1.,
        #                             dispersion_maj_b = disn_maj_b,
        #                             amplitude_min_b = 1.,
        #                             dispersion_min_b = disn_min_b,
        #                             signum_maj_b = sign_maj_b,
        #                             signum_min_b = sign_min_b, pa_b =
        #                             pa_b, dtype = target.dtype,
        #                             centering = 'origin', forreal =
        #                             True)
        #
        ikernel = self._igaussian_2dp(naxis1 = target.shape[1], naxis2
                                      = target.shape[0], cdelt1 = 1.,
                                      cdelt2 = 1., amplitude_maj_a =
                                      amp_maj_a, dispersion_maj_a =
                                      dis_maj_a, signum_maj_a = sign_maj_a,
                                      amplitude_min_a = amp_min_a,
                                      dispersion_min_a = dis_min_a,
                                      signum_min_a = sign_min_a, pa_a = pa_a,
                                      amplitude_maj_b = a_maj_b,
                                      dispersion_maj_b = dis_maj_b,
                                      amplitude_min_b = a_min_b,
                                      dispersion_min_b = dis_min_b,
                                      signum_maj_b = sign_maj_b,
                                      signum_min_b = sign_min_b, pa_b =
                                      pa_b, dtype = target.dtype,
                                      centering = 'origin', forreal =
                                      True)

        fft = pyfftw.builders.rfft2(target.copy(), planner_effort=None, threads= threads, auto_align_input=True, auto_contiguous=True, avoid_copy=False, norm=None)
        itarget = fft()

        iconvolved = ikernel*itarget

        ifft = pyfftw.builders.irfft2(iconvolved, s=target.shape, planner_effort=None, threads=7, auto_align_input=True, auto_contiguous=True, avoid_copy=False, norm=None)
        convolved = ifft()
        
        hdu[0].data[:] = convolved.astype(hdu[0].data.dtype)

        xmin = np.cos(pa_b)*HPBW_min_b/2
        ymin = np.sin(pa_b)*HPBW_min_b/2
        xmaj = -np.sin(pa_b)*HPBW_maj_b/2
        ymaj = np.cos(pa_b)*HPBW_maj_b/2
        
        print('Final result at centre', hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2])
        print('Final result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)])
        print('Final result at half power', hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)])
        hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2] += 2
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)] += 2.
        hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)] += 2.
        hdu.writeto(reconvolve_output_image, overwrite = True)
        print('Image to be found at', reconvolve_output_image)
        hdu.close()
        print('')
        return
    
    def _reconvolve(self, originplane, targetplane, originbeam, targetbeam, threads = 1):
        """(De-)convolve an image from an original to a targt resolution

        Input:
        originplane (ndarray): Original plane
        targetplane (ndarray): Target plane
        originbeam  (ndarray): Beam specifications (major axis
                               dispersion, minor axis dispersion,
                               position angle in rad, sin position
                               angle, cos position angle) assumed for
                               original
        targetbeam  (ndarray): Beam specifications (major axis
                               dispersion, minor axis dispersion,
                               position angle in rad, sin position
                               angle, cos position angle) assumed for
                               target
        threads (int)        : Numbers of CPUs to use
        
        
        Determines the Fourier transform FT(Go) of the Gaussian kernel
        Go that is described by the triplet originbeam and the Fourier
        transform FT(Gt) of the Gaussian kernel Gt described by
        targetbeam. The FFT of originplane is then divided by FT(G1)
        and multiplied by FT(G2) and transformed back into the image
        plane. A dispersion of np.inf in this context means a
        constant, a dispersion of 0 a delta function (equivalently a
        constant in the Fourier domain).

        """
        ikernel = self._igaussian_2dp(naxis1 = targetplane.shape[1],
                                      naxis2 = targetplane.shape[0],
                                      dispersion_maj_a = originbeam[0],
                                      signum_maj_a = 1.,
                                      dispersion_min_a = originbeam[1],
                                      signum_min_a = 1.,
                                      pa_a = originbeam[2],
                                      dispersion_maj_b = targetbeam[0],
                                      signum_maj_b = -1.,
                                      dispersion_min_b = targetbeam[1],
                                      signum_min_b = -1.,
                                      pa_b = targetbeam[2],
                                      dtype = targetplane.dtype,
                                      centering = 'origin',
                                      forreal = True)

        # Notice the use of threads as being passed, not from the instance
        fft = pyfftw.builders.rfft2(originplane, planner_effort =
                                    None, threads = threads,
                                    auto_align_input = True,
                                    auto_contiguous = True, avoid_copy
                                    = False, norm = None)
        itargetplane = fft()

        iconvolved = ikernel*itargetplane

        # Notice the use of threads as being passed, not from the instance
        ifft = pyfftw.builders.irfft2(iconvolved, s =
                                      targetplane.shape,
                                      planner_effort = None, threads =
                                      threads, auto_align_input =
                                      True, auto_contiguous = True,
                                      avoid_copy = False, norm = None)
        targetplane[:] = ifft()
        return

    def createstcubes(self, gauprops = [], outcubi = [], naxis = 4, naxis1 = 257, naxis2 = 513, pixelsize = 8.3333335E-04, ctype3 = 'VRAD', channelwidth = 5000, cellscal = None, restfreq = None, bmaj = None, bmin = None, bpa = None, noputinheader = [], overwrite = True):
        """Create a set of cubes to test beach

        Input:
        gauprops (list of ndarrays)   : Properties of Gaussians generated
        outcubi (list of str)         : Names of output cubes
        naxis  (int)                  : Number of axes
        naxis1 (int)                  : Size of cubes axis 1
        naxis2 (int)                  : Size of cubes axis 2
        pixelsize (float)             : Spacial pixel size in deg
        ctype3 (string)               : Type of third axis ('VRAD', 'FREQ')
        channelwidth (float)          : Channel width in m/s
        cellscal (None type or str)   : CELLSCAL, if None, no CELLSCAL
        restfreq (None type or float) : Rest frequency
        bmaj (None type or float)     : BMAJ, if None, no BMAJ
        bmin (None type or float)     : BMIN, if None, no BMIN
        bpa (None type or float)      : BPA, if None, no BPA
        noputinheader (list of lists) : List of lists of elements that
                                        should not appear in the
                                        headers
        overwrite (bool)              : Overwrite existing cubes?
                                        True: yes

        gauprops is a list of properties of single plane images of the
        cubes generated. The number of list elements is the number of
        cubes produced. Their names are listed in outcubi (which has
        to have as many elements as gauprops). The method produces
        cubes with one Gaussian in each plane. Each element of
        gauprops is an 2D ndarray with size (planes, 6), where planes
        is the number of planes, and the columns denote in that order,
        central position x in pixels (starting at 0), central position
        in y in pixels (starting at 0), amplitude, beam major HPBW in
        pixels, beam minor HPBW in pixels, beam position angle. Each
        cube has the same number naxis of axes, RA, DEC, VRAD (if
        naxis > 2), STOKES (if naxis > 3). If naxis < 3 (we hence
        produce not a cube but only one image), the elements of
        gauprops should accordingly have only one row. The data type
        is always intensity with the unit Jy/beam, the projection type
        is always sine projection J2000, and the cubes are centered on
        RA = 60 deg and DEC = -30 deg. The numbers of pixels in RA and
        DEC direction is the same for all cubes, in ctype3 direction
        it is determined by the number of rows in the corresponding
        elements in gauprops, the length of the STOKES axis is always
        1. The pixel size in RA and DEC direction is always the same
        and determined by the parameter pixelsize. The type of the
        third axis is determined by the keyword ctype3. The channel
        width is the same for all cubes and determined by the
        parameter channelwidth. The keyword cellscal ('CONSTANT' or
        '1/F') can be set by the parameter cellscal (if set to None,
        no keyword will appear in the header). Parameters restfreq
        (rest frequency in Hz), bmaj (average beam major axis HPBW in
        deg), bmin (average beam minor axis HPBW in deg), bpa (average
        beam position angle in deg) set the corresponding keywords in
        the headers. Choose '1420405745.51' for restfreq if you want
        the HI rest frequency.

        For each Gaussian generated with gauprops, the keywords BMAJi,
        BMINi, BPAi are generated in the header with the measures as
        specified in gauprops, unless they appear in the noputinheader
        list. noputinheader is a list of lists of strings. Each list
        in noputinheader specifies the keywords (works only for
        'BMAJi', 'BMINi', 'BPAi') that should not appear in the
        corresponding header.
        
        Finally, if overwrite is set to True,
        then the output cubes will be replaced if they exist,
        otherwise an error is thrown.

        """
        if len(gauprops) == 0:
            return
        if len(gauprops) != len(outcubi):
            return

        for i in range(len(gauprops)):
            if naxis == 4:
                newar = np.zeros((1, gauprops[i].shape[0], naxis2, naxis1), dtype = '>f4')
            elif naxis == 3 or gauprops[i].shape[0] != 1:
                naxis = 3
                newar = np.zeros((gauprops[i].shape[0], naxis2, naxis1), dtype = '>f4')
            elif naxis == 2:
                newar = np.zeros((naxis2, naxis1), dtype = '>f4')

            if naxis > 2:
                newar[newar.shape[0]//2,newar.shape[1]//2] = 1.
            else:
                newar[newar.shape[0]//2,newar.shape[1]//2] = 1.

            # Create a hdulist
            hdu = fits.HDUList([fits.PrimaryHDU(newar)])
            target = hdu[0].data.astype('float'+'{:d}'.format(hdu[0].data.itemsize*8))
            hdu[0].header['BSCALE'] = 1
            hdu[0].header['BZERO'] = 0
            hdu[0].header['BUNIT'] = 'JY/BEAM'
            hdu[0].header['BTYPE'] = 'intensity'
            hdu[0].header['CRPIX1'] = naxis1//2+1
            hdu[0].header['CDELT1'] = -pixelsize
            hdu[0].header['CRVAL1'] = 60.
            hdu[0].header['CTYPE1'] = 'RA---SIN'
            hdu[0].header['CRPIX2'] = naxis1//2+1
            hdu[0].header['CDELT2'] = pixelsize
            hdu[0].header['CRVAL2'] = -30.
            hdu[0].header['CTYPE2'] = 'DEC--SIN'
            hdu[0].header['OBJECT'] = 'p{:3d}'.format(i)
            hdu[0].header['EPOCH'] = '2000'
            hdu[0].header['EQUINOX'] = 2000.
            hdu[0].header['ORIGIN'] = 'beach'
            if naxis > 2:
                hdu[0].header['CRPIX3'] = 1
                hdu[0].header['CDELT3'] = channelwidth
            if 'FREQ' in ctype3:
                if type(restfreq) != type(None):
                    hdu[0].header['CRVAL3'] = restfreq
                else:
                    hdu[0].header['CRVAL3'] = self.HIFREQ                    
            hdu[0].header['CTYPE3'] = ctype3
            if naxis > 3:
                hdu[0].header['CRPIX4'] = 1
                hdu[0].header['CDELT4'] = 1
                hdu[0].header['CRVAL4'] = 1
                hdu[0].header['CTYPE4'] = 'STOKES'
            if type(cellscal) != type(None):
                if cellscal == 'constant':
                    hdu[0].header['CELLSCAL'] = 'CONSTANT'
                else:
                    hdu[0].header['CELLSCAL'] = '1/F'
            if type(restfreq) != type(None):
                hdu[0].header['RESTFREQ'] = restfreq
            if type(bmaj) != type(None):
                hdu[0].header['BMAJ'] = bmaj
            if type(bmaj) != type(None):
                hdu[0].header['BMIN'] = bmin
            if type(bpa) != type(None):
                hdu[0].header['BPA'] = bpa

            # Generate Gaussians and put their properties in header
            for j in range(gauprops[i].shape[0]):
                darray = self._gaussian_2dp( naxis1 = naxis1,
                                             naxis2 = naxis2, cdelt1 =
                                             1., cdelt2 = 1.,
                                             amplitude_maj_a =
                                             gauprops[i][j,2],
                                             dispersion_maj_a =
                                             gauprops[i][j,3]/np.sqrt(np.log(256.)),
                                             signum_maj_a = -1,
                                             amplitude_min_a = 1.,
                                             dispersion_min_a =
                                             gauprops[i][j,4]/np.sqrt(np.log(256.)),
                                             signum_min_a = -1, pa_a =
                                             np.pi*gauprops[i][j,5]/180.,
                                             dtype = target.dtype,
                                             centering =
                                             [gauprops[i][j,0],gauprops[i][j,1]],
                                             forreal = False).astype(hdu[0].data.dtype)
                if naxis == 2:
                    hdu[0].data[:] = darray.astype(hdu[0].data.dtype)
                if naxis == 3:
                    hdu[0].data[j,:] = darray.astype(hdu[0].data.dtype)
                if naxis == 4:
                    hdu[0].data[0,j,:] = darray.astype(hdu[0].data.dtype)


                if naxis > 2:
                    hdu[0].header['BMAJ{:d}'.format(j+1)] = gauprops[i][j,3]*pixelsize
                    hdu[0].header['BMIN{:d}'.format(j+1)] = gauprops[i][j,4]*pixelsize
                    hdu[0].header['BPA{:d}'.format(j+1)] = gauprops[i][j,5]
                else:
                    hdu[0].header['BMAJ'] = gauprops[i][j,3]*pixelsize
                    hdu[0].header['BMIN'] = gauprops[i][j,4]*pixelsize
                    hdu[0].header['BPA'] = gauprops[i][j,5]
                        
            
            hdu.writeto(outcubi[i], overwrite = overwrite)

def testplot():
    # Log-Normal Distribution

    mu, sigma = 0, 0.5

    measured = np.random.lognormal(mu, sigma, 1000)
    hist, edges = np.histogram(measured, density=True, bins=50)
    
    x = np.linspace(0.0001, 8.0, 1000)
    pdf = 1/(x* sigma * np.sqrt(2*np.pi)) * np.exp(-(np.log(x)-mu)**2 / (2*sigma**2))
    #cdf = (1+special.erf((np.log(x)-mu)/(np.sqrt(2)*sigma)))/2

    r = []
    for i in range(4):
        r += [ bokeh_plotting.figure(title='Log Normal Distribution', tools='pan,box_zoom,wheel_zoom,reset,save',
                               background_fill_color='#fafafa') ]
        r[-1].plot_height = 300
        r[-1].plot_width = 400
        r[-1].quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
           fill_color='navy', line_color='white', alpha=0.5)
        r[-1].line(x, pdf, line_color='#ff8888', line_width=4, alpha=0.7, legend_label='PDF')
        r[-1].y_range.start = 0
        r[-1].legend.location = 'center_right'

    q = bokeh_plotting.gridplot([[r[0], r[1]], [r[2], r[3]]])
    bokeh_plotting.output_file('Test.html')
    bokeh_plotting.save(q)

    # There is no straightforward method to export plots in a simple
    # way with bokeh, therefore one needs to add a matplotlib version
    # as well

       
def printcubeinfo(cubename):
    print()
    print('###################################')
    print('Header info',cubename,':')
    print('###################################')
    print()
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
    print()
    print('###############################')
    print(' Input and output structures :')
    print('###############################')
    print()

    print('Input cube 1')
    print('bmaj      : ', beach.binfo_input[0][:,0])
    print('bmin      : ', beach.binfo_input[0][:,1])
    print('bpa       : ', beach.binfo_input[0][:,2])
    print('freq      : ', beach.binfo_input[0][:,3])
    print('dx        : ', beach.binfo_input[0][:,4])
    print('rfreq/freq: ', beach.binfo_input[0][:,5])
    print('bmaj*freq : ', beach.binfo_input[0][:,6])
    print('bmin*freq : ', beach.binfo_input[0][:,7])
    print('Input cube 2')
    print('bmaj      : ', beach.binfo_input[1][:,0])
    print('bmin      : ', beach.binfo_input[1][:,1])
    print('bpa       : ', beach.binfo_input[1][:,2])
    print('freq      : ', beach.binfo_input[1][:,3])
    print('dx        : ', beach.binfo_input[1][:,4])
    print('rfreq/freq: ', beach.binfo_input[1][:,5])
    print('bmaj*freq : ', beach.binfo_input[1][:,6])
    print('bmin*freq : ', beach.binfo_input[1][:,7])
    print()
    print('Pixel cube 1')
    print('bmaj      : ', beach.binfo_pixel[0][:,0])
    print('bmin      : ', beach.binfo_pixel[0][:,1])
    print('bpa       : ', beach.binfo_pixel[0][:,2])
    print('freq      : ', beach.binfo_pixel[0][:,3])
    print('dx        : ', beach.binfo_pixel[0][:,4])
    print('rfreq/freq: ', beach.binfo_pixel[0][:,5])
    print('bmaj*freq : ', beach.binfo_pixel[0][:,6])
    print('bmin*freq : ', beach.binfo_pixel[0][:,7])
    print('Pixel cube 2')
    print('bmaj      : ', beach.binfo_pixel[1][:,0])
    print('bmin      : ', beach.binfo_pixel[1][:,1])
    print('bpa       : ', beach.binfo_pixel[1][:,2])
    print('freq      : ', beach.binfo_pixel[1][:,3])
    print('dx        : ', beach.binfo_pixel[1][:,4])
    print('rfreq/freq: ', beach.binfo_pixel[1][:,5])
    print('bmaj*freq : ', beach.binfo_pixel[1][:,6])
    print('bmin*freq : ', beach.binfo_pixel[1][:,7])
    print()
    print('Target cube 1')
    print('bmaj      : ', beach.binfo_target[0][:,0])
    print('bmin      : ', beach.binfo_target[0][:,1])
    print('bpa       : ', beach.binfo_target[0][:,2])
    print('freq      : ', beach.binfo_target[0][:,3])
    print('dx        : ', beach.binfo_target[0][:,4])
    print('rfreq/freq: ', beach.binfo_target[0][:,5])
    print('bmaj*freq : ', beach.binfo_target[0][:,6])
    print('bmin*freq : ', beach.binfo_target[0][:,7])
    print('Target cube 2')
    print('bmaj      : ', beach.binfo_target[1][:,0])
    print('bmin      : ', beach.binfo_target[1][:,1])
    print('bpa       : ', beach.binfo_target[1][:,2])
    print('freq      : ', beach.binfo_target[1][:,3])
    print('dx        : ', beach.binfo_target[1][:,4])
    print('rfreq/freq: ', beach.binfo_target[1][:,5])
    print('bmaj*freq : ', beach.binfo_target[1][:,6])
    print('bmin*freq : ', beach.binfo_target[1][:,7])
    print()

if __name__ == '__main__':

    #Beach().convoltests()

    print('')
    print('##################')
    print('##################')
    print('##################')
    print(' Final Tests/Demo ')
    print('##################')

    # The following parameters will be re-used in many tests
    naxis1 = 257
    naxis2 = 513
    naxis3 = 4

    # The Gaussian properties increase by one per plane and the
    # position angle by 5 per plane
    amp0 = 1.
    ainc = 0.0
    pinc = 5.
    bmaj0 = 10. 
    bmin0 = 7.
    binc = 0.
    bpa0 = 0
    cinc = 0
    
    gauprops = []
    
    # 2 cubes
    for i in range(2):
        planar = np.zeros((naxis3,6))
        for j in range(naxis3):
            planar[j, 0] = naxis1//2+(i*naxis3+j)*pinc
            planar[j, 1] = naxis2//2+(i*naxis3+j)*pinc
            planar[j, 2] = amp0+(i*naxis3+j)*ainc
            planar[j, 3] = bmaj0+(i*naxis3+j)*binc
            planar[j, 4] = bmin0+(i*naxis3+j)*binc
            planar[j, 5] = bpa0+(i*naxis3+j)*cinc
        gauprops.append(planar)
        
    print()
    print('########################')
    print('########################')
    print(' Test 1: Frequency axis')
    print('########################')
    print()

    incubi = ['test1_incubus1.fits', 'test1_incubus2.fits']
    Beach(verb = False).createstcubes(gauprops = gauprops, outcubi = incubi, naxis1 = naxis1, naxis2 = naxis2, ctype3='FREQ')
    print('Created input cubes {:s} and {:s}'.format(incubi[0], incubi[1]))
    printcubeinfo(incubi[0])
    printcubeinfo(incubi[1])
    outcubi = ['test1_outcubus1.fits', 'test1_outcubus2.fits']

    # This is just a preliminary exercise to create a parameter input interface
    params = { 'cubes': incubi,
               'tra_fitsnames': outcubi
    }
    beach = Beach(cubenames = params['cubes'], gentrans_exe = False)
    printbeachconts(beach)
    print()
    print('########################')
    
    beach.gentrans(tra_fitsnames = params['tra_fitsnames'], tra_overwrite = True)
    print()
    print('Created output cubes {:s} and {:s}'.format(outcubi[0], outcubi[1]))
    print()
    print('########################')
    print('########################')
    print(' Test 2: Velocity axis')
    print('########################')
    print()

    incubi = ['test2_incubus1.fits', 'test2_incubus2.fits']
    Beach(verb = False).createstcubes(gauprops = gauprops, outcubi = incubi, naxis1 = naxis1, naxis2 = naxis2, ctype3='VRAD')
    print('Created input cubes {:s} and {:s}'.format(incubi[0], incubi[1]))
    printcubeinfo(incubi[0])
    printcubeinfo(incubi[1])
    outcubi = ['test2_outcubus1.fits', 'test2_outcubus2.fits']

    # This is just a preliminary exercise to create a parameter input interface
    params = { 'cubes': incubi,
               'tra_fitsnames': outcubi
    }
    beach = Beach(cubenames = params['cubes'], gentrans_exe = False)
    printbeachconts(beach)
    beach.gentrans(tra_fitsnames = params['tra_fitsnames'], tra_overwrite = True)
    print()
    print('Created output cubes {:s} and {:s}'.format(outcubi[0], outcubi[1]))
    print()
    print('########################')
    print('########################')
    print('Test 3: Short version (shorter is still possible)')
    print('########################')
    print()

    incubi =  ['test2_incubus1.fits', 'test2_incubus2.fits']
    outcubi = ['test3_outcubus1.fits', 'test3_outcubus2.fits']
    params = { 'cubes': incubi,
               'tra_fitsnames': outcubi
    }
    Beach(cubenames = params['cubes'], tra_fitsnames = params['tra_fitsnames'], tra_overwrite = True)
    print()
    print('Created output cubes {:s} and {:s}'.format(outcubi[0], outcubi[1]))
    print()
    print('########################')
    print('########################')
    print(' Test 4: Just images')
    print('########################')
    print()
    gauprops[0] = gauprops[0][0,:].reshape((1,6))
    gauprops[1] = gauprops[1][0,:].reshape((1,6))
    incubi = ['test4_incubus1.fits', 'test4_incubus2.fits']
    Beach(verb = False).createstcubes(gauprops = gauprops, outcubi = incubi, naxis = 2, naxis1 = naxis1, naxis2 = naxis2, ctype3='VRAD')
    print('Created input cubes {:s} and {:s}'.format(incubi[0], incubi[1]))
    printcubeinfo(incubi[0])
    printcubeinfo(incubi[1])
    outcubi = ['test4_outcubus1.fits', 'test4_outcubus2.fits']

    # This is just a preliminary exercise to create a parameter input interface
    params = { 'cubes': incubi,
               'tra_fitsnames': outcubi
    }
    beach = Beach(cubenames = params['cubes'], gentrans_exe = False)
    printbeachconts(beach)
    beach.gentrans(tra_fitsnames = params['tra_fitsnames'], tra_overwrite = True)
    print()
    print('Created output cubes {:s} and {:s}'.format(outcubi[0], outcubi[1]))
    print()

#    testplot()
    