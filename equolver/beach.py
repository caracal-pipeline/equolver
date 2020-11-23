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
import datetime
import bokeh.plotting as bokeh_plotting
import bokeh.models as bokeh_models
import bokeh.layouts as bokeh_layouts
import bokeh.io as bokeh_io
import argparse
import sys
import textwrap

version = '0.0.4'

class Beach:
    """
    Lame acronym for a utility to equalize synthesized beams over
    channels

    Class variables:
    HIFREQ: rest frequency of the HI line
    
    """
    HIFREQ = 1420405751.7667

    def __init__(self, inc_cubes = None,
                       bin_bmaj = np.nan, bin_bmaj_replace = False,
                       bin_bmin = np.nan, bin_bmin_replace = False,
                       bin_bpa = np.nan, bin_bpa_replace = False,
                       bin_restfreq = HIFREQ, bin_restfreq_replace = False,
                       bin_normfreq = 1.4E9,
                       bst_parameter = 'all',
                       bst_scaling = 'all',
                       bst_stype = 'all',
                       bst_sample = 'all', bst_percents = 90,
                       bst_tolerance = 0.001, bst_nsamps = 200, bst_epsilon = 0.0005,
                       bst_maxiter = 1000000,
                       hist_plotname = None, hist_interactive = None, 
                       hist_sample = 'total', hist_scaling = 'frequency',
                       hist_n_per_bin = 5, hist_overwrite = False,
                       tar_bmaj_inter = ['bmaj', 'frequency', 'percentile', 'total'],
                       tar_bmaj_slope = 0.0,
                       tar_bmaj_absc = 0.0,
                       tar_bmin_inter = ['bmaj', 'frequency',  'percentile', 'total'],
                       tar_bmin_slope = 0.0,
                       tar_bmin_absc = 0.0,
                       tar_bpa_inter =  ['bpa', 'frequency',  'average', 'total'],
                       tar_bpa_slope =  0.0,
                       tar_bpa_absc =  0.0,
                       tar_scaling = 'frequency', genbstats_exe = True,
                       gentarget_exe = True,
                       gentrans_exe = True, tra_modelnames = None, tra_fitsnames = None,
                       tra_mode = 'mask', tra_hdmode = True, tra_tol = 0.0001,
                       tra_overwrite = False, tra_maxker = 1E7,
                       tra_commonbeam = True, tra_indibeam = True, tra_return_astropy = False,
                       threads = 1, verb = False):
        """Private instance variables:
        (multiple: None, a float, a list of floats, a numpy array, 
         or a list with numpy arrays)

        _inc_cubes (list)     : Names of input cubes
        _headers (list)         : List of open header hdus

        _bin_bmaj     (multiple)  : Beam major axis default value(s)
        _bin_bmaj_replace (bool)  : Enforce usage of default values?
                                  (True = yes)
        _bin_bmin     (multiple)  : Beam minor axis default value(s)
        _bin_bmin_replace (bool)  : Enforce usage of default values?
                              (True = yes)
        _bpa     (multiple)   : Beam position angle default value(s)
        _bpa_replace (bool)   : Enforce usage of default values? 

        _bin_restfreq (multiple)  : Rest frequency default value(s)
        _bpa_replace (bool)   : Enforce usage of default values? 

        _binfo_input (list) : List of arrays of point spread func-
                              tion, in the order of headernames and
                              headers each array of size (chans, 8),
                              where
 
                              first column bmaj in deg, 
                              
                              second column bmin in deg
                              
                              third column bpa in deg
                              
                              fourth column frequency in Hz

                              fifth column pixel size in deg, 

                              sixth column reference frequency divided
                              by frequency
          
                              seventh column bmaj in deg
                              scaled with frequency divided by the
                              normalisation frequency
 
                              eighth column bmin in deg scaled with
                              frequency divided by the normalisation
                              frequency

                              ninth column beam solid angle
                              (2pi/ln(256))*first_column*second_columnn,

                              10th column beam solid angle at
                              normalisation frequency
                              (2pi/ln(256))*seventh_column*eighth_columnn

                              eleventh column HPBW of a circular
                              Gaussian with BSA of column 9

                              twelfth column HPBW of a circular
                              Gaussian with BSA of column 10. chans is
                              the number of channels _binfo_pixel
                              (list) : _binfo_input converted into
                              pixel scaling using dispersion instead
                              of HPBW _bstats (dict) : Dictionary
                              containing all statistics _binfo_target
                              (list) : Target beam properties

        """
        self._initvars()
        self._verb = verb
        self._threads = threads
        self.genincubus(inc_cubes = inc_cubes)
        #self._initbinfo_inputvar(bin_bmaj = bin_bmaj, bin_bmaj_replace =
        #                         bin_bmaj_replace, bin_bmin = bin_bmin,
        #                         bin_bmin_replace = bin_bmin_replace,
        #                         bin_bpa = bin_bpa, bin_bpa_replace =
        #                         bin_bpa_replace, bin_restfreq =
        #                         bin_restfreq, bin_restfreq_replace =
        #                         bin_restfreq_replace, verb =
        #                         self._verb)
        
        for para in ['bin_bmaj', 'bin_bmaj_replace', 'bin_bmin',
                     'bin_bmin_replace', 'bin_bpa',
                     'bin_bpa_replace', 'bin_restfreq',
                     'bin_restfreq_replace']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])

        if genbstats_exe:
            self._genbstats_exe = True

        self.genbinfo(verb = self._verb)
        self._bin_normfreq = copy.deepcopy(bin_normfreq)

        for para in ['bst_parameter', 'bst_scaling', 'bst_stype', 'bst_sample',
                     'bst_percents', 'bst_tolerance', 'bst_nsamps', 'bst_epsilon',
                     'bst_maxiter', 'gentarget_exe']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])

        if self._genbstats_exe:
            self.genbstats(verb = self._verb)

        for para in ['hist_plotname', 'hist_interactive',
                     'hist_overwrite', 'hist_sample',
                     'hist_scaling', 'hist_overwrite', 'hist_n_per_bin']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])

        self.genhistoplots(verb = self._verb)

        for para in [ 'tar_bmaj_inter', 'tar_bmaj_slope',
                      'tar_bmaj_absc', 'tar_bmin_inter',
                      'tar_bmin_slope', 'tar_bmin_absc',
                      'tar_bpa_inter', 'tar_bpa_slope', 'tar_bpa_absc',
                      'tar_scaling','gentrans_exe']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])
            
        if self._gentarget_exe:
            self.gentarget(verb = self._verb)
            
        for para in [ 'gentrans_exe', 'tra_modelnames',
                      'tra_fitsnames', 'tra_mode', 'tra_hdmode', 'tra_tol',
                      'tra_maxker', 'tra_overwrite', 'tra_commonbeam',
                      'tra_indibeam', 'tra_return_astropy']:
            self.__dict__['_'+para] = copy.deepcopy(locals()[para])
            
        if self._gentrans_exe:
            self.gentrans(verb = self._verb)
            
        return

    def _initvars(self):
        """
        Reset/init instance variables
        """

        # Headers
        self._inc_cubes = None
        self._headers = None

        # Defaults
        self._bin_bmaj = None
        self._bin_bmaj_replace = False
        self._bin_bmin = None
        self._bin_bmin_replace = False
        self._bin_bpa = None
        self._bin_bpa_replace = False
        self._bin_restfreq = None
        self._bin_restfreq_replace = False

        # This is the normalisation for frequency-dependent beam
        # properties
        self._bin_normfreq = None

        # bmaj, bmin, pa, nu, pixelsize
        self._binfo_input = None

        # bmaj, bmin, sin pa, cos pa, nu
        self._binfo_pixel = None

        self._bst_parameter = None
        self._bst_scaling = None
        self._bst_stype = None
        self._bst_sample = None
        self._bst_percents = None
        self._bst_tolerance = None
        self._bst_nsamps = None
        self._bst_epsilon = None
        self._bst_maxiter = None

        # Statistics
        self._bstats = None

        self._hist_plotname = None
        self._hist_interactive = None
        self._hist_overwrite = None
        self._hist_sample = None
        self._hist_scaling = None
        self._hist_overwrite = None
        self._hist_n_per_bin = None
        
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

        self._tra_modelnames = None
        self._tra_fitsnames = None
        self._tra_mode = None
        self._tra_hdmode = None
        self._tra_tol = None
        self._tra_maxker = None
        self._tra_overwrite = None
        self._tra_commonbeam = None
        self._tra_indibeam = None
        self._tra_return_astropy = None

        self.tra_astropy = None

        self._threads = None
        self._verb = True
        return

    @property
    def inc_cubes(self):
        """
        Return a copy of inc_cubes
        """
        if type(self._inc_cubes) == type(None):
            return None
        return copy.deepcopy(self._binfo_input)

    @inc_cubes.setter
    def inc_cubes(self, value):

        # There is no use in letting the user change the cube names
        # but not changing the cubes, so this is enforced
        self.genincubus(inc_cubes = value)
        return
    
    @inc_cubes.deleter
    def inc_cubes(self, value):
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
    def bin_bmaj(self):
        """
        Return a copy of bin_bmaj
        """
        return copy.deepcopy(self._bin_bmaj)

    @bin_bmaj.setter
    def bin_bmaj(self, value):
        """
        Set bin_bmaj
        """
        self._bin_bmaj = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @bin_bmaj.deleter
    def bin_bmaj(self, value):
        self._bin_bmaj = np.nan
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def bin_bmaj_replace(self):
        return self._bin_bmaj_replace

    @bin_bmaj_replace.setter
    def bin_bmaj_replace(self, value):
        self._bin_bmaj_replace = value
    
    @bin_bmaj_replace.deleter
    def bin_bmaj_replace(self, value):
        self._bin_bmaj_replace = False

    @property
    def bin_bmin(self):
        """
        Return a copy of bin_bmin
        """
        return copy.deepcopy(self._bin_bmin)

    @bin_bmin.setter
    def bin_bmin(self, value):
        """
        Set bin_bmin
        """
        self._bin_bmin = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @bin_bmin.deleter
    def bin_bmin(self, value):
        self._bin_bmin = np.nan
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def bin_bmin_replace(self):
        return self._bin_bmin_replace

    @bin_bmin_replace.setter
    def bin_bmin_replace(self, value):
        self._bin_bmin_replace = value
    
    @bin_bmin_replace.deleter
    def bin_bmin_replace(self, value):
        self._bin_bmin_replace = False

    @property
    def bin_bpa(self):
        """
        Return a copy of bin_bpa
        """
        return copy.deepcopy(self._bin_bpa)

    @bin_bpa.setter
    def bin_bpa(self, value):
        """
        Set bin_bpa
        """
        self._bin_bpa = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @bin_bpa.deleter
    def bin_bpa(self, value):
        self._bin_bpa = np.nan
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def bin_bpa_replace(self):
        return self._bin_bpa_replace

    @bin_bpa_replace.setter
    def bin_bpa_replace(self, value):
        self._bin_bpa_replace = value
    
    @bin_bpa_replace.deleter
    def bin_bpa_replace(self, value):
        self._bin_bpa_replace = False

    @property
    def bin_restfreq(self):
        """
        Return a copy of bin_restfreq
        """
        return copy.deepcopy(self._bin_restfreq)

    @bin_restfreq.setter
    def bin_restfreq(self, value):
        """
        Set bin_restfreq
        """
        self._bin_restfreq = copy.deepcopy(value)
        self.genbinfo(verb = False)
        return

    @bin_restfreq.deleter
    def bin_restfreq(self, value):
        self._bin_restfreq = None
        self._binfo_input = None
        self._binfo_pixel = None
        self._bstats = None
        self._binfo_target = None
        return
        
    @property
    def bin_restfreq_replace(self):
        return self._bin_restfreq_replace

    @bin_restfreq_replace.setter
    def bin_restfreq_replace(self, value):
        self._bin_restfreq_replace = value
    
    @bin_restfreq_replace.deleter
    def bin_restfreq_replace(self, value):
        self._bin_restfreq_replace = False

    @property
    def bin_normfreq(self):
        """
        Return a copy of bin_normfreq
        """
        return copy.deepcopy(self._bin_normfreq)

    @bin_normfreq.setter
    def bin_normfreq(self, value):
        """
        Set bin_normfreq
        """
        self._bin_normfreq = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @bin_normfreq.deleter
    def bin_normfreq(self, value):
        self._bin_normfreq = 1.4E9
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def bst_parameter(self):
        """
        Return a copy of bst_parameter
        """
        return copy.deepcopy(self._bst_parameter)

    @bst_parameter.setter
    def bst_parameter(self, value):
        """
        Set bst_parameter
        """
        self._bst_parameter = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_parameter.deleter
    def bst_parameter(self, value):
        self._bst_parameter = 'all'
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def bst_scaling(self):
        """
        Return a copy of bst_scaling
        """
        return copy.deepcopy(self._bst_scaling)

    @bst_scaling.setter
    def bst_scaling(self, value):
        """
        Set bst_scaling
        """
        self._bst_scaling = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_scaling.deleter
    def bst_scaling(self, value):
        self._bst_scaling = 'frequency'
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def bst_stype(self):
        """
        Return a copy of bst_stype
        """
        return copy.deepcopy(self._bst_stype)

    @bst_stype.setter
    def bst_stype(self, value):
        """
        Set bst_stype
        """
        self._bst_stype = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_stype.deleter
    def bst_stype(self, value):
        self._bst_stype = ['median','medstdev']
        self._bstats = None
        self._binfo_target = None
        return
    
    @property
    def bst_sample(self):
        """
        Return a copy of bst_sample
        """
        return copy.deepcopy(self._bst_sample)

    @bst_sample.setter
    def bst_sample(self, value):
        """
        Set bst_sample
        """
        self._bst_sample = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_sample.deleter
    def bst_sample(self, value):
        self._bst_sample = 'total'
        self._bstats = None
        self._binfo_target = None
        return
    
    @property
    def bst_percents(self):
        """
        Return a copy of bst_percents
        """
        return copy.deepcopy(self._bst_percents)

    @bst_percents.setter
    def bst_percents(self, value):
        """
        Set bst_percents
        """
        self._bst_percents = copy.deepcopy(value)
        if self_genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_percents.deleter
    def bst_percents(self, value):
        self._bst_percents = 90
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def bst_tolerance(self):
        """
        Return a copy of bst_tolerance
        """
        return copy.deepcopy(self._bst_tolerance)

    @bst_tolerance.setter
    def bst_tolerance(self, value):
        """
        Set bst_tolerance
        """
        self._bst_tolerance = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_tolerance.deleter
    def bst_tolerance(self, value):
        self._bst_tolerance = 0.01
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def bst_nsamps(self):
        """
        Return a copy of bst_nsamps
        """
        return copy.deepcopy(self._bst_nsamps)

    @bst_nsamps.setter
    def bst_nsamps(self, value):
        """
        Set bst_nsamps
        """
        self._bst_nsamps = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_nsamps.deleter
    def bst_nsamps(self, value):
        self._bst_nsamps = 200
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def bst_epsilon(self):
        """
        Return a copy of bst_epsilon
        """
        return copy.deepcopy(self._bst_epsilon)

    @bst_epsilon.setter
    def bst_epsilon(self, value):
        """
        Set bst_epsilon
        """
        self._bst_epsilon = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_epsilon.deleter
    def bst_epsilon(self, value):
        self._bst_epsilon = 0.0005
        self._bstats = None
        self._binfo_target = None
        return

    @property
    def bst_maxiter(self):
        """
        Return a copy of bst_maxiter
        """
        return copy.deepcopy(self._bst_maxiter)

    @bst_maxiter.setter
    def bst_maxiter(self, value):
        """
        Set bst_maxiter
        """
        self._bst_maxiter = copy.deepcopy(value)
        if self._genbstats_exe:
            self.genbstats(verb = False)
        return

    @bst_maxiter.deleter
    def bst_maxiter(self, value):
        self._bst_maxiter = 1000000
        self._bstats = None
        self._binfo_target = None
        return

    ####
    
    @property
    def hist_interactive(self):
        """
        Return a copy of hist_interactive
        """
        return copy.deepcopy(self._hist_interactive)

    @hist_interactive.setter
    def hist_interactive(self, value):
        """
        Set hist_interactive
        """
        self._hist_interactive = copy.deepcopy(value)
        #self.genhistoplots(verb = False)
        return

    @hist_interactive.deleter
    def hist_interactive(self, value):
        self._hist_interactive = None
        return
    
    @property
    def hist_n_per_bin(self):
        """
        Return a copy of hist_n_per_bin
        """
        return copy.deepcopy(self._hist_n_per_bin)

    @hist_n_per_bin.setter
    def hist_n_per_bin(self, value):
        """
        Set hist_n_per_bin
        """
        self._hist_n_per_bin = copy.deepcopy(value)
        #self.genhistoplots(verb = False)
        return

    @hist_n_per_bin.deleter
    def hist_n_per_bin(self, value):
        self._hist_n_per_bin = 5
        return
    
    @property
    def hist_plotname(self):
        """
        Return a copy of hist_plotname
        """
        return copy.deepcopy(self._hist_plotname)

    @hist_plotname.setter
    def hist_plotname(self, value):
        """
        Set hist_plotname
        """
        self._hist_plotname = copy.deepcopy(value)
        #self.genhistoplots(verb = False)
        return

    @hist_plotname.deleter
    def hist_plotname(self, value):
        self._hist_plotname = None
        return
    
    @property
    def hist_sample(self):
        """
        Return a copy of hist_sample
        """
        return copy.deepcopy(self._hist_sample)

    @hist_sample.setter
    def hist_sample(self, value):
        """
        Set hist_sample
        """
        self._hist_sample = copy.deepcopy(value)
        #self.genhistoplots(verb = False)
        return

    @hist_sample.deleter
    def hist_sample(self, value):
        self._hist_sample = 'total'
        return
    
    @property
    def hist_scaling(self):
        """
        Return a copy of hist_scaling
        """
        return copy.deepcopy(self._hist_scaling)

    @hist_scaling.setter
    def hist_scaling(self, value):
        """
        Set hist_scaling
        """
        self._hist_scaling = copy.deepcopy(value)
        self.genhistoplots(verb = False)
        return

    @hist_interactive.deleter
    def hist_scaling(self, value):
        self._hist_scaling = 'frequency'
        return

    @property
    def hist_overwrite(self):
        """
        Return a copy of hist_overwrite
        """
        return copy.deepcopy(self._hist_overwrite)

    @hist_overwrite.setter
    def hist_overwrite(self, value):
        """
        Set hist_overwrite
        """
        self._hist_overwrite = copy.deepcopy(value)
        #self.genhistoplots(verb = False)
        return

    @hist_interactive.deleter
    def hist_overwrite(self, value):
        self._hist_overwrite = False
        return
###

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

    @tra_fitsnames.setter
    def tra_fitsnames(self, value):
        """
        Set tra_fitsnames
        """
        self._tra_fitsnames = copy.deepcopy(value)
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_fitsnames.deleter
    def tra_fitsnames(self):
        self._tra_fitsnames = None
        return

    @property
    def tra_modelnames(self):
        """
        Return a copy of _tra_modelnames
        """
        return self._tra_modelnames

    @tra_modelnames.setter
    def tra_modelnames(self, value):
        """
        Set tra_modelnames
        """
        self._tra_modelnames = copy.deepcopy(value)
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_modelnames.deleter
    def tra_modelnames(self):
        self._tra_modelnames = None
        return

    @property
    def tra_mode(self):
        """
        Return a copy of _tra_mode
        """
        return self._tra_mode

    @tra_mode.setter
    def tra_mode(self, value):
        """
        Set tra_mode
        """
        self._tra_mode = copy.deepcopy(value)
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_mode.deleter
    def tra_mode(self):
        self._tra_mode = 'mask'
        return

    @property
    def tra_hdmode(self):
        """
        Return a copy of _tra_hdmode
        """
        return self._tra_hdmode

    @tra_hdmode.setter
    def tra_hdmode(self, value):
        """
        Set tra_hdmode
        """
        self._tra_hdmode = copy.deepcopy(value)
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_hdmode.deleter
    def tra_hdmode(self):
        self._tra_hdmode = True
        return

    @property
    def tra_tol(self):
        """
        Return a copy of _tra_tol
        """
        return self._tra_tol

    @tra_tol.setter
    def tra_tol(self, value):
        """
        Set tra_tol
        """
        self._tra_tol = copy.deepcopy(value)
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_tol.deleter
    def tra_tol(self):
        self._tra_tol = 2.
        return
    
    @property
    def tra_maxker(self):
        """
        Return a copy of _tra_maxker
        """
        return self._tra_maxker

    @tra_maxker.setter
    def tra_maxker(self, value):
        """
        Set tra_maxker
        """
        self._tra_maxker = copy.deepcopy(value)
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_maxker.deleter
    def tra_maxker(self):
        self._tra_maxker = 1E7
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
    def tra_return_astropy(self):
        """
        Return a copy of tra_return_astropy
        """
        return self._tra_return_astropy

    @tra_return_astropy.setter
    def tra_return_astropy(self, value):
        """
        Set tra_return_astropy
        """
        self._tra_return_astropy = value
        if self._gentrans_exe:
            self.gentrans(verb = False)
        return

    @tra_return_astropy.deleter
    def tra_return_astropy(self):
        self._tra_return_astropy = False
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
        self._inc_cubes = None
        return

    def genincubus(self, inc_cubes = None):
        """
        Add inc_cubes to intrinsic inc_cubes

        Input:
        inc_cubes (list of str): List of names of input cubes/images

        Reads inc_cubes in as the list of target headers/images
        """
        if type(inc_cubes) == type(None):
            return

        # Check if inc_cubes have been defined before
        #dontinit = True
        #for i in range(len(inc_cubes)):
        #    if inc_cubes[i] != self._inc_cubes[i]:
        #        dontinit == False
        #        break
            
        # Save the time to reload the cubes    
        #if dontinit and silent:
        #    return
        
        self._resetheaders()
        self._inc_cubes = copy.deepcopy(inc_cubes)
        self._headers = []

        if type(self._inc_cubes) == type([]):
            cubenamelist = self._inc_cubes
        else:
            cubenamelist = [self._inc_cubes]    
            
        for cube in self._inc_cubes:
            if type(cube) == type(''):
                self._headers += [fits.getheader(cube)]
            else:
                self._headers += [cube[0].header]
            

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
            self.genincubus()
            
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
        

            self._binfo_input.append(np.empty((naxis3,12,)))
            self._binfo_input[-1][:] = np.nan
        return

    def _getdefault(self, inquant, quantname = None):
        """Fill target quantname, wich is the name of an instance
        variable, with the content of inquant or return formatted
        inquant
        
        Input:
        (multiple: None, a float, a list of floats, a numpy array, 
         or a list with numpy arrays)
        inquant (multiple): input
        quantname (str)   : Name of quantity (without leading _)

        The output is assumed to be a list of linear ndarrays of the
        size of naxis3 of self._headers (i.e. the number of channels),
        which should be filled with the default values for bin_bmaj,
        bmin, etc., following the same expansion scheme. If inquant is
        None, the call is ignored, if it is a single (float) number
        (including np.nan), all fields in the output are filled with
        that number, if it is a list of numbers (including np.nan)
        each ndarray in the output list is filled with the
        corresponding number, if it is an ndarray, every array in the
        output list is filled with this ndarray, and finally, if a
        list of ndarrays is provided, the output will completely be
        filled with the input. An input list has always to have as
        many elements as self._headers and an ndarray always has to
        have the same dimension as the naxis3 of the corresponding
        cube. If quantname is not none, the instance variable
        self._quantname is filled with the output.

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

    def _setdefreplace(self, bin_bmaj_replace = None,
                            bin_bmin_replace = None,
                            bin_bpa_replace = None,
                            bin_restfreq_replace = None):
        """
        Read in values for replacements of defvalues
        """
        
        if self._bin_bmaj_replace == None:
            self._bin_bmaj_replace = False
        if self._bin_bmin_replace == None:
            self._bin_bmin_replace = False
        if self._bin_bpa_replace == None:
            self._bin_bpa_replace = False
        if self._bin_restfreq_replace == None:
            self._bin_restfreq_replace = False

        if bin_bmaj_replace != None:
            self._bin_bmaj_replace = bin_bmaj_replace
        if bin_bmin_replace != None:
            self._bin_bmin_replace = bin_bmin_replace
        if bin_bpa_replace != None:
            self._bin_bpa_replace = bin_bpa_replace
        if bin_restfreq_replace != None:
            self._bin_restfreq_replace = bin_restfreq_replace
            
    def _initbinfo_inputvar(self, bin_bmaj = None, bin_bmaj_replace = None, bin_bmin = None,
                    bin_bmin_replace = None, bin_bpa = None, bin_bpa_replace = None,
                            bin_restfreq = None, bin_restfreq_replace = None,
                            bin_normfreq = None, verb = False):
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
    
    def genbinfo(self, bin_bmaj = None, bin_bmaj_replace = None,
                       bin_bmin = None, bin_bmin_replace = None,
                       bin_bpa = None,  bin_bpa_replace = None,
                       bin_restfreq = HIFREQ, bin_restfreq_replace = None,
                       bin_normfreq = None, verb = True):
        """
        Fill beam properties into info table
        
        Input:
        (multiple: None, a float, a list of floats, a numpy
                 array, or a list with numpy arrays)

        bin_bmaj     (multiple): Beam major axis default value(s)
        bin_bmaj_replace (bool): Enforce usage of default values?
                                 (True = yes)
        bin_bmin     (multiple): Beam minor axis default value(s)
        bin_bmin_replace (bool): Enforce usage of default values?
                                 (True = yes)
        bin_bpa     (multiple) : Beam position angle default value(s)
        bin_bpa_replace (bool) : Enforce usage of default values? 

        bin_restfreq (multiple)    : Rest frequency default value(s)
        bin_restfreq_replace (bool): Enforce usage of default values?

        bin_normfreq (float)       : Frequency to normalize beam to if
                                 mode is 'frequency'

        Note that if None is passed as a value, the input is ignored
        if self._quantity is not None.
        """
        stop = self._initbinfo_inputvar(bin_bmaj = bin_bmaj,
                                        bin_bmaj_replace = bin_bmaj_replace,
                                        bin_bmin = bin_bmin, bin_bmin_replace =
                                        bin_bmin_replace, bin_bpa = bin_bpa,
                                        bin_bpa_replace = bin_bpa_replace,
                                        bin_restfreq = bin_restfreq,
                                        bin_restfreq_replace =
                                        bin_restfreq_replace,
                                        bin_normfreq = bin_normfreq, verb = verb)

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
                warning('binfo_input not present. Use genincubus first.')
            return
        
        # Do not try to regenerate, if this is none, it has deliberately been
        # set to None, as __init__ has a default of 1E9
        if type(self._bin_normfreq) == type(None):
            if verb or self._verb:
                warnings.warn('No normalization frequency read in, which '+ \
                              'disables further processing.')
            return

        print('genbinfo: reading beam info.')
        print()
        
        normfreq = self._getdefault(self._bin_normfreq)

        # Get defaults
        bmaj = self._getdefault(self._bin_bmaj)
        bmin = self._getdefault(self._bin_bmin)
        bpa = self._getdefault(self._bin_bpa)
        restfreq = self._getdefault(self._bin_restfreq)

        #self._setdefreplace(bin_bmaj_replace = bin_bmaj_replace,
                            #bmin_replace = bmin_replace,
                            #bpa_replace = bpa_replace,
                            #restfreq_replace = restfreq_replace)

        for i in range(len(self._headers)):

            # It would be rather surprising if there is a channel-dependent
            # rest frequency, but who knows
            restfreqtab = self._getchanval(
                'RESTFREQ', self._headers[i], value = restfreq[i],
                usevalue = bin_restfreq_replace, usedefault = True)
            
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
                 usevalue = bin_bmaj_replace, usedefault = True, dscal = self._binfo_input[i][:,5])
            self._binfo_input[i][:,1] = self._getchanval( \
                'BMIN', self._headers[i], value = bmin[i],
                usevalue = bin_bmin_replace, usedefault = True, dscal = self._binfo_input[i][:,5])
            self._binfo_input[i][:,2] = self._getchanval( \
                'BPA', self._headers[i], value = bpa[i],
                usevalue = bin_bpa_replace, usedefault = True)

            if np.isnan(self._binfo_input[i][:,3].sum()):
                thafreq = normfreq[i]
            else:
                thafreq = self._binfo_input[i][:,3]
            self._binfo_input[i][:,6] = self._binfo_input[i][:,0]*\
                thafreq/normfreq[i]
            self._binfo_input[i][:,7] = self._binfo_input[i][:,1]*\
                thafreq/normfreq[i]
            self._binfo_input[i][:,8] = (2*np.pi/np.log(256))*\
                self._binfo_input[i][:,0]*self._binfo_input[i][:,1]
            self._binfo_input[i][:,9] = (2*np.pi/np.log(256))*\
                self._binfo_input[i][:,6]*self._binfo_input[i][:,7]
            self._binfo_input[i][:,10] = np.sqrt(self._binfo_input[i][:,0]*self._binfo_input[i][:,1])
            self._binfo_input[i][:,11] = np.sqrt(self._binfo_input[i][:,6]*self._binfo_input[i][:,7])
            
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

            # bsa simply scale by pixel size
            boutfarray[:,8] = binfarray[:,8]/(binfarray[:,4]*binfarray[:,4])
            boutfarray[:,9] = binfarray[:,9]/(binfarray[:,4]*binfarray[:,4])

            # average beam
            boutfarray[:,10] = binfarray[:,10]/binfarray[:,4]/np.sqrt(np.log(256))
            boutfarray[:,11] = binfarray[:,11]/binfarray[:,4]/np.sqrt(np.log(256))
            
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
        for key0 in ['bmaj', 'bmin', 'bpa', 'bsa', 'ceb']:
            self._bstats[key0] = {}
            
            # scaling
            for key1 in ['constant', 'frequency']:
                self._bstats[key0][key1] = {}

                # bst_stype
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
    
    def _getcollist(self, sca, sam):
        """
        Helper function to genbstats or genhistoplots
        """
        
        # collect sample as list of np arrays
        collist = []

        
        if sam == 'cube':
            for i in self._binfo_input:
                if sca == 'constant':
                    collist.append(i[:,[0,1,2,8,10]])
                else:
                    collist.append(i[:,[6,7,2,9,11]])

        elif sam == 'chan':
            for i in range(self._bstats['bmaj']['constant']['maximum'][sam].shape[0]):
                collist.append(np.empty((0,5,), dtype=np.float64))
                for j in self._binfo_input:
                    if j.shape[0] > i:
                        if sca == 'constant':
                            collist[i] = np.append(collist[i],
                                               j[i:i+1,[0,1,2,8,10]],
                                               axis = 0)
                        else:
                            collist[i] = np.append(collist[i],
                                               j[i:i+1,[6,7,2,9,11]],
                                               axis = 0)
        elif sam == 'total':
            collist = [np.empty((0,5,))]
            for i in self._binfo_input:
                if sca == 'constant':
                    collist[0] = np.append(collist[0], i[:,[0,1,2,8,10]],
                                           axis = 0)
                else:
                    collist[0] = np.append(collist[0], i[:,[6,7,2,9,11]],
                                           axis = 0)
        else:
            return None            
        
        return collist
    
    def _initbstatsvar(self, bst_parameter = None, bst_stype = None,
                       bst_scaling = None, bst_sample = None,
                       bst_percents = None, bst_tolerance = None,
                       bst_nsamps = None, bst_epsilon = None,
                       bst_maxiter = None, verb = False):
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
            
    def genbstats(self, bst_parameter = None, bst_scaling = None,
                  bst_stype = None, bst_sample = None, bst_percents =
                  None, bst_tolerance = None, bst_nsamps = None,
                  bst_epsilon = None, bst_maxiter = None, verb = True):
        """
        Generate statistics and dump it into the bstats structure
        
        Input:
        bst_parameter (str or list of str): Parameter name ('bmaj', 
                                            'bmin', 'bpa', 'bsa', 'ceb')
        bst_scaling (str or list of str)  : Scaling type ('constant', 
                                        'frequency')
        bst_stype (str or list of str)    : Type of statistics to
                                            calculate ('minimum',
                                            'maximum', 'average',
                                            'stdev', 'median', 'mad',
                                            'madstdev', 'percentile',
                                            'percents', 'commonbeam')
        bst_sample (str or list of str)   : Sample(s) to calculate
                                            statistics on ('cube', 'chan',
                                            'total')
        bst_percents (float)              : Percents for the percentile
                                            statistics
        bst_tolerance (float)                 : Tolerance for the common beam
        bst_nsamps (int)                      : Number of edges of beam for
                                            common beam
        bst_epsilon (float)                   : Epsilon for common beam
        bst_maxiter (int)                   : Maximum iterations for common beam

        The method generates statistics on the collected beam
        properties.  The parameters bst_parameter, bst_scaling, bst_sample, and
        type determine which part of the bstats structure gets
        filled. If for any of the parameters 'all' is chosen (which is
        the default), all fields are filled. If the scaling type is
        'constant', the given values for bmaj and bmin are evaluated,
        if 'frequency' is chosen, all beam sizes are scaled to the
        same norm-frequency nf, assuming that the beam sizes scale
        with 1/frequency. b(nf) = b(f)*f/nf . bsa is the beam solid
        angle of a beam, ceb the circular equivalent beam, the
        circular beam with a beam solid angle of bsa (sqrt(bmaj*bmin))

        The _bstats entity is a nested dictionary
        _bstats[bst_parameter][bst_scaling][bst_stype][bst_sample], where sample
        determines the type of element (numpy scalar, numpy ndarray,
        or list of scalars, see below):

        bst_parameter: 
            major axis dispersions/hpbws ('bmaj')
            minor axis dispersions/hpbws ('bmin')
            beam position angles ('bpa')
            beam solid angle ('bsa')
            circular equivalent beam dispersions/hpbws('ceb')
        bst_scaling:
            constant ('const')
            frequency ('frequency')
        bst_stype:
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
            Parameters bst_tolerance, bst_nsamps, epsilon are used for this 
            method ('commonbeam')

        bst_sample:
            statistics to be carried out for all channels per cube
            ('cube', generates lists with length of the number of
            input cubes)

            statistics to be carried out for all cubes per channel
            ('chan', generates linear ndarrays with a length of the
            maximum number of channels in any cube)

            statistics to be carried out for all channels in all
            cubes ('total', generates a float)

        """

        stop = self._initbstatsvar(bst_parameter = bst_parameter,
                                   bst_scaling = bst_scaling,
                                   bst_stype = bst_stype, bst_sample =
                                   bst_sample, bst_percents =
                                   bst_percents, bst_tolerance =
                                   bst_tolerance, bst_nsamps =
                                   bst_nsamps, bst_epsilon =
                                   bst_epsilon, bst_maxiter =
                                   bst_maxiter, verb = verb)

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

        print('genbstats: generating statistics')
        
        para = copy.deepcopy(self._bst_parameter)
        scal = copy.deepcopy(self._bst_scaling)
        styp = copy.deepcopy(self._bst_stype)
        samp = copy.deepcopy(self._bst_sample)

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
            para = ['bmaj', 'bmin', 'bpa', 'bsa', 'ceb']
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
                warnings.warn('Input information struct '+ \
                              'cannot be generated. Returning without '+ \
                              'generating statistics.')
            return

        # This is needed later, map for parameters
        parma = {'bmaj': 0, 'bmin': 1, 'bpa': 2, 'bsa': 3, 'ceb': 4}

        # Now let's do it

        # get struct
        #struct = self._binfo_input
        for sca in scal:
            
            for sam in samp:
                
                collist = self._getcollist(sca, sam)

                if type(collist) == type(None):
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
                        'percents' : self._bst_percents,
                        'tolerance': self._bst_tolerance,
                        'nsamps'   : self._bst_nsamps,
                        'epsilon'  : self._bst_epsilon,
                        'maxiter'  : self._bst_maxiter
                        }
                    stats = met(collist, **kwargs)
                    #print('stats ', stats)

                    # Get function to generate statistics
                    for par in para:

                        #print('par', par, 'parma', parma[par])

                        # Put the right value into array
                        self._bstats[par][sca][sty][sam] = stats[:,parma[par]]
                        #print('self._bstats[', par, '][', sty, '][', sam,'] =', self._bstats[par][sca][sty][sam])

                        if sty != 'percents':
                            print('Statistics Parameter: {}\n'.format(par)+\
                                  'Scaling:              {}\n'.format(sca)+\
                                  'Statistics:           {}\n'.format(sty)+\
                                  'Sample:               {}'.format(sam))
                            if sam == 'cube':
                                for i in range(len(self._bstats[par][sca][sty][sam])):
                                    print('  cube {}: {}'.format(i,self._bstats[par][sca][sty][sam][i]))
                            if sam == 'chan':
                                for i in range(self._bstats[par][sca][sty][sam].shape[0]):
                                    print('  channel {}: {}'.format(i,self._bstats[par][sca][sty][sam][i]))
                            if sam == 'total':
                                print('  total: {}'.format(self._bstats[par][sca][sty][sam][0]))
                            print()
        print()
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
        
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(output, np.amin(array, axis = 0).reshape((1,5,)),
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
        
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(output, np.amax(array, axis = 0).reshape((1,5,)),
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
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(output, np.mean(array, axis = 0).reshape((1,5,)),
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

        output = np.empty((0, 5))
        for array in parray:
            
            # Make sure that there is no error if there is only one
            if array.shape[0] < 2:
                dof = 0.
            else:
                dof = 1.
        
            output = np.append(output, np.std(array, axis = 0,
                                      ddof = dof).reshape((1,5,)), axis = 0)

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
        
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(output, np.median(array, axis = 0).reshape((1,5,)),
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
        
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(output, stats.median_abs_deviation(
                array, axis = 0, scale = 1.).reshape((1,5,)), axis = 0)

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
        
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(output, stats.median_abs_deviation(
                array, axis = 0, scale = 'normal').reshape((1,5,)), axis = 0)

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
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(output, np.percentile(
                array, kwargs['percents'], axis = 0).reshape((1,5,)), axis = 0)

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
        
        output = np.empty((0, 5))
        for array in parray:
            output = np.append(
                output, np.array([kwargs['percents']]).repeat(5).reshape((1,5,)),
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
        
        output = np.empty((0, 5))
        for array in parray:

            # Unit of this one does not matter
            bmaj=array[:,0] * u.deg

            # Unit of this one does not matter
            bmin=array[:,1] * u.deg

            # Unit of this one has to be deg
            bpa=array[:,2] * u.deg

            print('bmaj', bmaj)
            
            my_beams = radio_beam.Beams(bmaj, bmin, bpa)

            repeat = True
            tolerance = kwargs['tolerance']
            maxiter = kwargs['maxiter']
            while repeat:
                try:
                    common_beam = my_beams.common_beam(
                        tolerance = tolerance,
                        nsamps = kwargs['nsamps'],
                        epsilon = kwargs['epsilon'],
                        maxiter = maxiter
                    )
                    repeat = False
                except radio_beam.utils.BeamError:
                    tolerance = tolerance*1.1
                    maxiter = int(maxiter*1.1)
                    
            bmaj = common_beam.major.to(u.degree).value
            bmin = common_beam.minor.to(u.degree).value
            bpa = common_beam.pa.to(u.degree).value
            bsa = (2*np.pi/np.log(256))*bmaj*bmin
            ceb = np.sqrt(bmaj*bmin)
            
            output = np.append(
                output, np.array([bmaj, bmin, bpa, bsa, ceb]).reshape((1,5,)),
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
                        if type(self._bst_parameter) == type(''):
                            if self._bst_parameter != item[0]:
                                self._bst_parameter = [self._bst_parameter]
                        if type(self._bst_scaling) == type(''):
                            if self._bst_scaling != item[1]:
                                self._bst_scaling = [self._bst_scaling]
                        if type(self._bst_stype) == type(''):
                            if self._bst_stype != item[2]:
                                self._bst_stype = [self._bst_stype]
                        if type(self._bst_sample) == type(''):
                            if self._bst_sample != item[2]:
                                self._bst_sample = [self._bst_sample]
                        if not item[0] in self._bst_parameter:
                            self._bst_parameter += [item[0]]
                        if not item[1] in self._bst_scaling:
                            self._bst_scaling += [item[0]]
                        if not item[0] in self._bst_stype:
                            self._bst_stype += [item[0]]
                        if not item[0] in self._bst_sample:
                            self._bst_sample += [item[0]]
                        self.genbstats(bst_parameter = item[0], bst_scaling = item[1], bst_stype =
                                       item[2], bst_sample = item[3], verb = verb)
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
    
    def _gaussian(self, x, cent, amp, sigma):
        """
        Gaussian
        Input:
        cent (float): centre
        amp (float) : amplitude
        sigma (float) : sigma
        Return:
        gaussian() Gaussian
        """
        return amp*np.exp(-0.5*np.power((x-cent)/sigma,2))


    def _inithistovar(self, hist_plotname = None, hist_interactive =
                      None, hist_sample = None, hist_scaling = None,
                      hist_overwrite = None, hist_n_per_bin = None,
                      verb = True):
        output = False
        paras = locals().copy()
        paras.pop('self')
        paras.pop('verb')
        
        for param in paras.keys():
            if type(paras[param]) != type(None):
                self.__dict__['_'+param] = copy.deepcopy(paras[param])

        # Check if we do any plot
        if type(self._hist_plotname) == type(None) and type(self._hist_interactive) == type(None):
            return output
            
        # For the rest we need only some parameters to be defined
        paras.pop('hist_plotname')
        paras.pop('hist_interactive')

        for param in paras.keys():
            if self.__dict__['_'+param] == None:
                if self._verb or verb:
                    warnings.warn('Parameter {} is not defined.'.format(param))
                output = True
        return output

    def genhistoplots(self, hist_plotname = None, hist_interactive =
                      None, hist_sample = None, hist_scaling = None,
                      hist_overwrite = None, hist_n_per_bin = None, verb
                      = True):
        """
        Generate a histogram
        
        Input:
        hist_plotname (str)   : Name of plot
        hist_interactive (str): Name of interactive plot
        hist_sample (str)     : Sample to plot 'cube', 'chan', or 'total'
        hist_scaling (str)    : Scaling to use ('frequency' or 'constant')
        hist_overwrite (bool) : Allow overwriting files produced before
        verb (bool)           : Use verbose mode (if) True or not

        Generates histograms of the beam properties read in. This will
        either produce static histograms or an interactive html
        file. With hist_sample the user can choose which statistics
        should be generated. 'cube' means that the histogram is
        generated for each cube, 'chan' means that it is generated for
        a specific channel in all cubes, 'total' (the default) means
        all channels in all cubes. hist_scaling decides if the
        statistics are plotted using the original beam properties or
        if they are first scaled to the norm frequency assuming that
        the beam scales proportionally to the inverse of the
        frequency.

        """
        arguments = locals().copy()
        arguments.pop('self')
        stop = self._inithistovar(**arguments)
        if stop:
            if verb or self._verb:
                warnings.warn('Parameters missing. Not generating histogram.')
            return

        # If no plots are requested, just return
        if type(self._hist_plotname) == type(None) and type(self._hist_interactive) == type(None):
            return

        print('genhistoplots: generating diagnostic plots.')
        print()

        plotname = self._hist_plotname
        intername = self._hist_interactive
        
        if type(self._hist_plotname) == type('') and self._hist_plotname[-4:] != '.png':
            plotname = self._hist_plotname+'.png'
        if type(self._hist_interactive) == type('') and self._hist_interactive[-5:] != '.html':
            intername = self._hist_interactive+'.html'

        # Get the sample
        collist = self._getcollist(self._hist_scaling, self._hist_sample)

        # Get as much statistics as possible
        stats = []
        fivepars = ['bmaj', 'bmin', 'bpa', 'bsa', 'ceb']            
        fivetitles = ['BMAJ', 'BMIN', 'BPA', 'BSA', 'CEB']
        fiveunits = ['arcsec', 'arcsec', 'deg', 'sqarcsec', 'arcsec']
        fivescaling = [3600., 3600., 1., 3600*3600, 3600]
        children = []
        children2 = []
        #print('Got here scaling, sample', self._hist_scaling, self._hist_sample)
        #print('Got here stats intern:', self._bstats['bmaj']['constant']['average']['chan'])
        for i in range(len(collist)):
            k = len('{:d}'.format(len(collist)-1))
            
            if self._hist_sample == 'total':
                titlestring = '<b>Beam properties of all planes in all data sets</b>'
            elif self._hist_sample == 'cube':
                if type(self._inc_cubes[i]) == type(''):
                    titlestring = '<b>Beam properties of planes in data set {:s}</b>'.format(self._inc_cubes[i])
                else:
                    titlestring = '<b>Beam properties of planes in data set {:d}</b>'.format(i)
                    
            elif self._hist_sample == 'chan':
                titlestring = '<b>Beam properties of channel {:d}</b>'.format(i)
            else:
                if verb or self._verb:
                    warnings.warn('hist_sample not properly defined. Not generating histogram.')
                return              

            bins = -(-collist[i][:,0].size//self._hist_n_per_bin)
            fivehists = [
                np.histogram(collist[i][:,0]*3600., bins = bins),
                np.histogram(collist[i][:,1]*3600., bins = bins),
                np.histogram(collist[i][:,2], bins = bins),
                np.histogram(collist[i][:,3]*3600.*3600., bins = bins),
                np.histogram(collist[i][:,4]*3600., bins = bins)                
            ]
            r = []
                        
            for j in range(5):
                
                # Get the stats if available
                r += [ bokeh_plotting.figure(tools=['pan','box_zoom','wheel_zoom','reset,save'],
                                             background_fill_color='#fafafa', x_axis_label='{:s} ({:s})'.format(fivetitles[j], fiveunits[j]), y_axis_label = 'Counts') ]
                r[-1].plot_height = 300
                r[-1].plot_width = 300
                r[-1].y_range.start = 0

                # Histogram
                r[-1].quad(top=fivehists[j][0], bottom=0, left=fivehists[j][1][:-1], right=fivehists[j][1][1:],
                           fill_color='navy', line_color='white', alpha=0.5)

                # Gaussian descriptors
                avgaux = np.linspace(fivehists[j][1][0], fivehists[j][1][-1], 200)
                amplitude = float(len(collist))/bins
                average = self._bstats[fivepars[j]][self._hist_scaling]['average'][self._hist_sample][i]*fivescaling[j]
                stdev = self._bstats[fivepars[j]][self._hist_scaling]['stdev'][self._hist_sample][i]*fivescaling[j]
                amplitudeave = collist[i][:,0].size*(fivehists[j][1][-1]-fivehists[j][1][0])/(np.sqrt(2.*np.pi)*stdev*bins)

                # Plot average and stdev properties if they exist
                if np.isfinite(average) and np.isfinite(stdev):
                    # Calculate Gaussian
                    avgauy = self._gaussian(avgaux, average, amplitudeave, stdev)

                    # Add Gaussian to the plot
                    gaussline = r[-1].line(avgaux, avgauy, line_color='#ff6666', line_width=4, alpha=0.7, legend_label='Av')
                    #r[-1].add_tools(bokeh_models.HoverTool(tooltips="""<font size-"3">Gaussian from<br>rms and stdev<br>$x $y</font>""", renderers=[gaussline], mode='vline'))
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips="""<font size-"3">Gaussian from<br>rms and stdev<br>$x $y</font>""", renderers=[gaussline]))

                    # Add lines representing average and stdev                       
                    gaussaverage = r[-1].line([average,average],[0,amplitudeave], line_color='#ff6666', line_width=4, alpha=0.7, legend_label='Av')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('average', '{:.1f}'.format(average))], attachment='right', renderers=[gaussaverage]))

                    # To avoid funny looking plots
                    avmst = average-stdev
                    if avgaux[0] > average-stdev:
                        avmst = avgaux[0]
                    avpst = average+stdev
                    if avgaux[-1] < average+stdev:
                        avpst = avgaux[-1]
                    gaussstdev = r[-1].line([avmst,avpst],[self._gaussian(average-stdev, average, amplitudeave, stdev),self._gaussian(average+stdev, average, amplitudeave, stdev)], line_color='#ff6666', line_width=4, alpha=0.7, legend_label='Av')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('stdev', '{:.1f}'.format(stdev))], attachment='right', renderers=[gaussstdev]))
                    
                # Gaussian from mad around median

                # Gaussian descriptors
                average = self._bstats[fivepars[j]][self._hist_scaling]['median'][self._hist_sample][i]*fivescaling[j]
                stdev = self._bstats[fivepars[j]][self._hist_scaling]['madstdev'][self._hist_sample][i]*fivescaling[j]
                amplitudemed = collist[i][:,0].size*(fivehists[j][1][-1]-fivehists[j][1][0])/(np.sqrt(2.*np.pi)*stdev*bins)

                # Plot average and stdev properties if they exist
                if np.isfinite(average) and np.isfinite(stdev):
                    # Calculate Gaussian
                    avgauy = self._gaussian(avgaux, average, amplitudemed, stdev)

                    # Add Gaussian to the plot
                    gaussline = r[-1].line(avgaux, avgauy, line_color='#66ffff', line_width=4, alpha=0.7, legend_label='Med')
                    #r[-1].add_tools(bokeh_models.HoverTool(tooltips="""<font size-"3">Gaussian from<br>rms and stdev<br>$x $y</font>""", renderers=[gaussline], mode='vline'))
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips="""<font size-"3">Gaussian from<br>rms and stdev<br>$x $y</font>""", renderers=[gaussline]))

                    # Add lines representing average and stdev                       
                    gaussaverage = r[-1].line([average,average],[0,amplitudemed], line_color='#66ffff', line_width=4, alpha=0.7, legend_label='Med')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('median', '{:.1f}'.format(average))], attachment='right', renderers=[gaussaverage]))

                    # To avoid funny looking plots
                    avmst = average-stdev
                    if avgaux[0] > average-stdev:
                        avmst = avgaux[0]
                    avpst = average+stdev
                    if avgaux[-1] < average+stdev:
                        avpst = avgaux[-1]
                    gaussstdev = r[-1].line([avmst,avpst],[self._gaussian(average-stdev, average, amplitudemed, stdev),self._gaussian(average+stdev, average, amplitudemed, stdev)], line_color='#66ffff', line_width=4, alpha=0.7, legend_label='Med')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('madstd', '{:.1f}'.format(stdev))], attachment='right', renderers=[gaussstdev]))

                # Vertical line dimensions (not using Spans because they don't have legends)
                #height = r[-1].__attr__
                #print('height',height)

                ymax = max(amplitudeave, amplitudemed, np.amax(fivehists[j][0]))
                
                # Min and max
                minimum = self._bstats[fivepars[j]][self._hist_scaling]['minimum'][self._hist_sample][i]*fivescaling[j]
                if np.isfinite(minimum):
                    miniline = r[-1].line([minimum,minimum],[0,ymax], line_color='#ff66b3', line_width=4, alpha=0.7, legend_label='Minmax')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('minimum', '{:.1f}'.format(minimum))], attachment='right', renderers=[miniline]))
                maximum = self._bstats[fivepars[j]][self._hist_scaling]['maximum'][self._hist_sample][i]*fivescaling[j]
                if np.isfinite(maximum):
                    maxiline = r[-1].line([maximum,maximum],[0,ymax], line_color='#ff66b3', line_width=4, alpha=0.7, legend_label='Minmax')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('maximum', '{:.1f}'.format(maximum))], attachment='left', renderers=[maxiline]))

                # Percentile
                perc = self._bstats[fivepars[j]][self._hist_scaling]['percentile'][self._hist_sample][i]*fivescaling[j]
                percent = self._bstats[fivepars[j]][self._hist_scaling]['percents'][self._hist_sample][i]
                if np.isfinite(perc) and np.isfinite(percent):
                    percline = r[-1].line([perc,perc],[0,ymax], line_color='#ffb366', line_width=4, alpha=0.7, legend_label='Perc')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('{}% percentile'.format(percent), '{:.1f}'.format(perc))], attachment='left', renderers=[percline]))
                
                # Add possibility to toggle graphs and lines
                r[-1].legend.click_policy="hide"

                # Common beam
                common = self._bstats[fivepars[j]][self._hist_scaling]['commonbeam'][self._hist_sample][i]*fivescaling[j]
                if np.isfinite(common):
                    commonline = r[-1].line([common,common],[0,ymax], line_color='#66ffb3', line_width=4, alpha=0.7, legend_label='Common')
                    r[-1].add_tools(bokeh_models.HoverTool(tooltips=[('common beam', '{:.1f}'.format(common))], attachment='right', renderers=[commonline]))

            children += [bokeh_models.widgets.markups.Div(text=titlestring, style={'font-size': '125%', 'font-family': 'bf', 'color': 'black'})]
            children += [bokeh_plotting.gridplot([[r[0], r[1], r[2], r[3], r[4]]])]
            
            children2 += [bokeh_models.widgets.markups.Div(text=titlestring, style={'font-size': '125%', 'font-family': 'bf', 'color': 'black'})]
            children2 += [bokeh_plotting.gridplot([[r[0], r[1], r[2], r[3], r[4]]], toolbar_location=None)]
            
        q = bokeh_layouts.column(children=children)
        q2 = bokeh_layouts.column(children=children2)
        if type(intername) != type(None):
            bokeh_plotting.output_file(intername)
            bokeh_plotting.save(q)
        if type(plotname) != type(None):
            bokeh_io.export_png(q2, filename=plotname)
    

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
        denoting (in that order) bst_parameter, bst_scaling,
        bst_stype, and bst_sample as described in method
        genbstats. The corresponding values will then be copied from
        the bstats struct.

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
        if type(self._bin_normfreq) == type(None):
            if verb or self._verb:
                warnings.warn('No normalization frequency read in, which '+ \
                              'disables further processing.')
            return

        normfreq = self._getdefault(self._bin_normfreq)

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

        print('gentarget: generating target beam properties')
        print()
        
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

    def _initgentransvar(self, tra_modelnames = None, tra_fitsnames =
                         None, tra_mode = None, tra_hdmode = None,
                         tra_tol = None, tra_maxker = None,
                         tra_commonbeam = None, tra_indibeam = None,
                         tra_overwrite = None, tra_return_astropy =
                         None, threads = None, verb = True):
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
        paras.pop('tra_modelnames')
        for param in paras.keys():
            if self.__dict__['_'+param] == None:
                if verb or self._verb:
                    warnings.warn('Parameter {} is not defined.'.format(param))
                output = True
        return output
            
    def gentrans(self, tra_modelnames = None, tra_fitsnames = None,
                 tra_mode = None, tra_hdmode = None, tra_tol = None,
                 tra_commonbeam = None, tra_indibeam = None,
                 tra_indimode = None, tra_overwrite = None, tra_return_astropy
                 = None, tra_maxker = None, threads = None, verb = True):
        """(De-)convolve input data cubes or images to target beam shapes

        Input:
        tra_modelnames (str or list of str): Input fits file names, 
                                             containing the models,
                                             alternatively astropy
                                             hdulists can be given.
                                             None is a valid input.
        tra_fitsnames (str or list of str) : Output fits file names
                                             None is a valid input.
        tra_mode (str)                     : 'scale', 'mask', 
                                             'hybrid', 'max'
        tra_tol (float)                    : tolerance to determine if
                                             convolution failed
        tra_commonbeam (bool)              : Generate common (average)
                                             beam information in 
                                             header
        tra_indibeam (bool)                : Generate individual beam 
                                             information in header
        tra_hdmode (bool)                  : Generate information about
                                             scaling/convolution in header
        tra_overwrite (bool)               : Overwrite output if
                                             already existent
                                             (True: yes)?
        tra_return_astropy (bool)                  : Return list of astropy hdulists?
                                             (True: yes)

        tra_maxker (float):                  Maximum value that the FT of the 
                                             convolution kernel can assume, will
                                             assume failure if larger
        threads (bool)                     : Number of threads

        Serially opens all cubes (images) listed in inputnames and
        generates cubes (De-)convolved to the resolution as listed in
        the target structure. It then convolves all cubes listed in
        tra_modelnames with the respective Gaussians and adds those to
        the output. The number of cubes and the dimensionality of the
        cubes must be identical to the one of the input cubes if not
        None. 

        A (de-)convolution is declared a success or a failure by
        comparing the sums of a plane divided by the beam. If the
        ratio of the larger sum divided by the smaller sum is larger
        than tra_tol, the (de-)convolution has failed. In addition,
        tra_maxker allows the user to manually set the maximum
        exponent in the FT of a kernel as an alternative way to decide
        whether a reconvolution is bound to fail.

        The mode of the deconvolution is determined by parameter
        tra_mode as follows:
          
          'scale': Do not convolve but scale the intensity to the
          target beam (divide by original beam solid angle and
          multiply with target beam solid angle.)

          'mask': (De-)convolve and mask channel in the output if the
          deconvolution fails. 

          'hybrid': Attempt to (de-) convolve the plane and fall back
          to scale if the (de-)convolution fails.

          'max': Attempt to (de-) convolve the plane. If that fails,
          convolve along the beam minor axis to the target minor beam
          if that is larger than the original. Then scale.

        If tra_hdmode is set to True, a keyword 'EQMODE' with the
        value of mode is added to the header. In addition, if EQMODE
        is not scale or mask, any plane for which the first
        convolution failed, is highlighted by the keyword-value pair
        EQS_i = 'SCALE' or EQSC = 'HYBRID'. i in this context is the
        plane number (Fortran/FITS style).

        """
        stop = self._initgentransvar(tra_modelnames = tra_modelnames,
                                     tra_fitsnames = tra_fitsnames,
                                     tra_mode = tra_mode, tra_tol =
                                     tra_tol, tra_overwrite =
                                     tra_overwrite, tra_commonbeam =
                                     tra_commonbeam, tra_indibeam =
                                     tra_indibeam, tra_return_astropy
                                     = tra_return_astropy, tra_maxker
                                     = tra_maxker, verb = verb,
                                     threads = threads)
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

        # Also check if tra_fitsnames have the same type and number as inputnames,
        # same for tra_modelnames
        
        if type(self._tra_fitsnames) == type(''):
            transn = [self._tra_fitsnames]
        else:
            transn = self._tra_fitsnames

        if type(self._inc_cubes) == type(''):
            cuben = [self._inc_cubes]
        else:
            cuben = self._inc_cubes

        if type(self._tra_modelnames) == type([]):
            modeln = self._tra_modelnames
        else:
            modeln = [self._tra_modelnames]

        if len(transn) != len(cuben):
            if verb or self._verb:
                warnings.warn('Number of input names not matching number of'+ \
                              'output data sets. Not generating output data.')
            return

        if type(modeln) != type(None): 
            if len(modeln) != len(cuben):
                if verb or self._verb:
                    warnings.warn('Number of input names not matching number of'+ \
                                  'output data sets. Not generating output data.')
                return
            
        if not self._tra_mode in ['scale', 'mask', 'hybrid', 'max']: 
            if verb or self._verb:
                warnings.warn('Number of input names not matching number of'+ \
                                  'output data sets. Not generating output data.')
            return
            

        # After this marginal verification, we continue
        print('gentrans: applying beams')
        print()

        # Check if the method is supposed to return a value
        if self._tra_return_astropy:
            tra_return_astropylist = []
        
        # Open, reconvolve, copy
        for i in range(len(cuben)):
            if type(cuben[i]) == type(''):
                incubus = fits.open(cuben[i])
            else:
                incubus = cuben[i]
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

            planeinfo = {}
            
            for plane in range(incubus_image.shape[0]):
                print()
                if type(cuben[i]) == type(''):
                    print('Processing {:s} plane {:d}'.format(cuben[i],plane))
                else:
                    print('Processing cube {:d} plane {:d}'.format(i,plane))                    
                originbeam = self._binfo_pixel[i][plane,:]
                targetbeam = self._binfo_target[i][plane,:]
                originplane = incubus_image[plane,:,:]
                targetplane = outcubus_image[plane,:,:]
                print('HPBW in pixels PA in deg, bmaj, bmin, bpa: {:.1f} {:.1f} {:.1f} -> {:.1f} {:.1f} {:.1f}'.format(originbeam[0]*np.sqrt(np.log(256.)),originbeam[1]*np.sqrt(np.log(256.)),180.*originbeam[2]/np.pi,targetbeam[0]*np.sqrt(np.log(256.)),targetbeam[1]*np.sqrt(np.log(256.)),180.*targetbeam[2]/np.pi))

                if self._tra_mode == 'scale':

                    # Scale by the area of the beam
                    targetplane[:] = originplane*targetbeam[0]*targetbeam[1]/(originbeam[0]*originbeam[1])

                else:
                    
                    # Attempt to convolve
                    self._reconvolve(originplane, targetplane,
                                     originbeam, targetbeam, threads =
                                     self._threads, maxker =
                                     self._tra_maxker)

                    failure = False
                    
                    # Assess success: one nan/inf pixel marks failure
                    if np.isfinite(targetplane).astype(int).sum() < targetplane.size:
                        failure = True
                    else:
                        # Asses success: the inner quarter of the image
                        # should approximately show the same sum,
                        # normalised by the beam size
                        #originflux = originplane[originplane.shape[0]//4:3*originplane.shape[0]//4,
                        #            originplane.shape[1]//4:3*originplane.shape[1]//4].sum()/(originbeam[0]*originbeam[1])
                        #targetflux = targetplane[targetplane.shape[0]//4:3*targetplane.shape[0]//4,
                        #            targetplane.shape[1]//4:3*targetplane.shape[1]//4].sum()/(targetbeam[0]*targetbeam[1])
                        originflux = originplane.sum()/(originbeam[0]*originbeam[1])
                        targetflux = targetplane.sum()/(targetbeam[0]*targetbeam[1])
                        failure = np.amax([np.fabs(originflux), np.fabs(targetflux)])/np.amin([np.fabs(originflux), np.fabs(targetflux)]) > self._tra_tol+1.

                    if failure:
                        if type(cuben[i]) == type(''):
                            print('{:s} plane {:d}: failed to re-convolve'.format(cuben[i],plane))
                        else:
                            print('Cube {:d} plane {:d}: failed to re-convolve'.format(i,plane))                            
                        
                        if self._tra_mode == 'mask':
                            print('Masking')
                            targetplane[:] = np.nan
                            
                        if self._tra_mode == 'hybrid':
                            print('Scaling instead')
                            targetplane[:] = originplane*targetbeam[0]*targetbeam[1]/(originbeam[0]*originbeam[1])
                            planeinfo['EQS_{:d}'.format(plane+1)] = 'SCALE'
                            #targetplane[:] = np.nan
                            
                        if self._tra_mode == 'max':
                            print('Max approach')

                            # Target beam is smaller or oriented differentely than origin beam
                            # Find out if bmin in targetbeam is larger than bmin in originbeam
                            if targetbeam[1] > originbeam[1]:

                                print('HPBW in pixels PA in deg, bmaj, bmin, bpa: {:.1f} {:.1f} {:.1f} -> {:.1f} {:.1f} {:.1f}'.format(originbeam[0]*np.sqrt(np.log(256.)),originbeam[1]*np.sqrt(np.log(256.)),180.*originbeam[2]/np.pi,originbeam[0]*np.sqrt(np.log(256.)),targetbeam[1]*np.sqrt(np.log(256.)),180.*originbeam[2]/np.pi))

                                # Then we can go half way, we convolve to the same minor beam,
                                # but not more
                                self._reconvolve(originplane,
                                                 targetplane,
                                                 originbeam,
                                                 [originbeam[0],
                                                  targetbeam[1],
                                                  originbeam[2]],
                                                 maxker =
                                                 self._tra_maxker,
                                                 threads =
                                                 self._threads)

                                # Check again
                                if np.isfinite(targetplane).astype(int).sum() < targetplane.size:
                                    print('This here')
                                    failure_again = True
                                else:
                                    # Asses success: the inner quarter of the image
                                    # should approximately show the same sum,
                                    # normalised by the beam size
                                    #orginflux = originplane[originplane.shape[0]//4:3* originplane.shape[0]//4,
                                    #            originplane.shape[1]//4:3*originplane.shape[1]//4].sum()/(originbeam[0]*originbeam[1])
                                    #targetflux = targetplane[targetplane.shape[0]//4:3*targetplane.shape[0]//4,
                                    #            targetplane.shape[1]//4:3*targetplane.shape[1]//4].sum()/(targetbeam[0]*targetbeam[1])
                                    originflux = originplane.sum()/(originbeam[1])
                                    targetflux = targetplane.sum()/(targetbeam[1])
                                    failure_again = np.amax([np.fabs(originflux), np.fabs(targetflux)])/np.amin([np.fabs(originflux), np.fabs(targetflux)]) > self._tra_tol+1.

                                if failure_again:
                                    if type(cuben[i]) == type(''):
                                        print('{:s} plane {:d}: failed to re-convolve again'.format(cuben[i],plane))
                                    else:
                                        print('Cube {:d} plane {:d}: failed to re-convolve again'.format(i,plane))
                                    print('Rescaling only.')
                                    
                                    # If we failed again (just a safeguard), we scale only
                                    targetplane[:] = originplane*targetbeam[0]*targetbeam[1]/(originbeam[0]*originbeam[1])
                                    planeinfo['EQS_{:d}'.format(plane+1)] = 'SCALE'

                                else:
                                    print('Applying hybrid approach')
                                    # If we succeeded we need to scale the rest
                                    targetplane[:] = targetplane[:]*targetbeam[0]/originbeam[0]
                                    planeinfo['EQS_{:d}'.format(plane+1)] = 'HYBRID'
                            else:
                                print('Reconvolving not possible, just scalingyo')
                                # If we succeeded we need to scale the rest
                                #targetplane[:] = targetplane*targetbeam[1]/originbeam[1]
                                targetplane[:] = originplane*targetbeam[0]*targetbeam[1]/(originbeam[0]*originbeam[1])
                                planeinfo['EQS_{:d}'.format(plane+1)] = 'HYBRID'                                

            if type(modeln[0]) != type(None):
                if type(modeln[i]) == type(''):
                    inmodel = fits.open(modeln[i])
                else:
                    inmodel = modeln[i]
                inmodel_image = inmodel[0].data.astype('float'+'{:d}'.format(inmodel[0].data.itemsize*8))
                modelshape = inmodel_image.shape
                if inmodel_image.ndim > 3:
                    inmodel_image = np.squeeze(inmodel_image)
                if inmodel_image.ndim == 2:
                    inmodel_image = inmodel_image.reshape((1, inmodel_image.shape[0], inmodel_image.shape[1]))
                outmodel_image = inmodel_image.copy()*0.+np.nan
                
                originbeam = [0.,0.,0.,0.,1.]
                for plane in range(incubus_image.shape[0]):
                    if type(cuben[i]) == type(''):
                        print('Processing {:s} plane {:d}'.format(modeln[i],plane))
                    else:
                        print('Processing model {:d} plane {:d}'.format(i,plane))                            
                    targetbeam = self._binfo_target[i][plane,:]
                    originplane = inmodel_image[plane,:,:]
                    targetplane = outmodel_image[plane,:,:]
                    self._reconvolve(originplane, targetplane, originbeam,
                                     targetbeam, threads = self._threads)
                    
                # Replace incubus with outcubus and write out
                incubus[0].data[:] = outcubus_image.astype(incubus[0].data.dtype)+\
                    outmodel_image.astype(incubus[0].data.dtype)
                
            else:
                # Replace incubus with outcubus and write out
                incubus[0].data[:] = outcubus_image.astype(incubus[0].data.dtype)

            if self._tra_hdmode:
                incubus[0].header['EQMODE'] = self._tra_mode
                for key in planeinfo.keys():
                    incubus[0].header[key] = planeinfo[key]
                
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

            # Document occurrence of equolver in history
            incubus[0].header['HISTORY'] = ''
            incubus[0].header['HISTORY'] = '{:s} Generated by equolver'.format(datetime.datetime.now().strftime("%d-%b-%Y %H:%M:%S"))
            incubus[0].header['HISTORY'] = 'See https://github.com/caracal-pipeline/equolver'
            incubus[0].header['HISTORY'] = ''

            # Finally write and close cube
            incubus.writeto(self._tra_fitsnames[i], overwrite = self._tra_overwrite)

            if self._tra_return_astropy:
                tra_return_astropylist += [incubus]
            else:
                incubus.close()

        if self._tra_return_astropy:
            self.tra_astropy = tra_return_astropylist
            return tra_return_astropylist
        else:
            return

    def _gaussian_2dp(self, naxis1 = 100, naxis2 = 100, cdelt1 = 1.,
                      cdelt2 = 1., amplitude_maj_a = 1.,
                      dispersion_maj_a = np.inf, signum_maj_a = -1.,
                      amplitude_min_a = 1., dispersion_min_a = np.inf,
                      signum_min_a = -1., pa_a = 0., amplitude_maj_b =
                      1., dispersion_maj_b = np.inf, amplitude_min_b =
                      1., dispersion_min_b = np.inf, signum_maj_b =
                      1., signum_min_b = 1., pa_b = 0., dtype =
                      'float32', centering = 'origin', forreal = True,
                      maxker = 1E300):
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
        maxker (float):           Maximum value in return array, returns 
                                  nan if any pixel gets larger than that.

        Calculates the product of four planar quasi-Gaussians

            P = amplitude_maj_a*amplitude_maj_b*
                exp(signum_maj_a/2*(x0a/dispersion_maj_a)+
                signum_min_a/2*(x1a/dispersion_min_a))*
                amplitude_min_a*amplitude_min_b*
                exp(signum_maj_b/2*(x0b/dispersion_maj_b)+
                signum_min_b/2*(x1b/dispersion_min_b))
                
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

        #exponent=exponent1a+exponent0a+exponent1b+exponent0b
        #if np.isnan(exponent.sum()):
        #    return exponent1a*np.nan
        #print(exponent.max())
        #if exponent.max() > maxex:
        #    return exponent1a*np.nan
        rval = amplitude*np.exp(0.5*(exponent1a+exponent0a+exponent1b+exponent0b))
        if rval.max() > maxker:
            return rval*np.nan
        return rval

    def _igaussian_2dp(self, naxis1 = 100, naxis2 = 100, cdelt1 = 1.,
                      cdelt2 = 1., amplitude_maj_a = 1.,
                      dispersion_maj_a = 0., signum_maj_a = -1.,
                      amplitude_min_a = 1., dispersion_min_a = 0.,
                      signum_min_a = -1., pa_a = 0., amplitude_maj_b =
                      1., dispersion_maj_b = 0., amplitude_min_b =
                      1., dispersion_min_b = 0., signum_maj_b =
                      1., signum_min_b = 1., pa_b = 0., dtype =
                       'float32', centering = 'origin', forreal = True,
                       maxker = 1E300):
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
                                  forreal, maxker = maxker)

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
    
    def _reconvolve(self, originplane, targetplane, originbeam,
                    targetbeam, maxker = 1E300, threads = 1):
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
        maxker (float):        Maximum value that the FT of the 
                               convolution kernel can assume, will
                               return a nan plane if larger
        
        
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
                                      forreal = True, maxker = maxker)

        if np.isnan(ikernel.sum()):
            targetplane[:] = np.nan
            return

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

    def createstcubes(self, mode = 'gauss', gauprops = [], outcubi = [], naxis = 4, naxis1 = 257, naxis2 = 513, pixelsize = 8.3333335E-04, ctype3 = 'VRAD', channelwidth = 5000, cellscal = None, restfreq = None, bmaj = None, bmin = None, bpa = None, noputinheader = [], overwrite = True):
        """Create a set of cubes to test beach

        Input:
        gauprops (list of ndarrays)   : Properties of Gaussians generated
        mode (str)                    : Produce Gaussians ('gauss') or point
                                        sources ('point')
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
        cubes with one Gaussian or one point source in each
        plane. Each element of gauprops is an 2D ndarray with size
        (planes, 6), where planes is the number of planes, and the
        columns denote in that order, central position x in pixels
        (starting at 0), central position in y in pixels (starting at
        0), amplitude, beam major HPBW in pixels, beam minor HPBW in
        pixels, beam position angle. If mode is 'point', then point
        sources are generated (but the beam properties put into the
        header), if it is 'gauss', then Gaussians are produced. Point
        sources are always generated exactly on top of a pixel using
        the nearest neighbour approach. Each cube has the same number
        naxis of axes, RA, DEC, VRAD (if naxis > 2), STOKES (if naxis
        > 3). If naxis < 3 (we hence produce not a cube but only one
        image), the elements of gauprops should accordingly have only
        one row. The data type is always intensity with the unit
        Jy/beam, the projection type is always sine projection J2000,
        and the cubes are centered on RA = 60 deg and DEC = -30
        deg. The numbers of pixels in RA and DEC direction is the same
        for all cubes, in ctype3 direction it is determined by the
        number of rows in the corresponding elements in gauprops, the
        length of the STOKES axis is always 1. The pixel size in RA
        and DEC direction is always the same and determined by the
        parameter pixelsize. The type of the third axis is determined
        by the keyword ctype3. The channel width is the same for all
        cubes and determined by the parameter channelwidth. The
        keyword cellscal ('CONSTANT' or '1/F') can be set by the
        parameter cellscal (if set to None, no keyword will appear in
        the header). Parameters restfreq (rest frequency in Hz), bmaj
        (average beam major axis HPBW in deg), bmin (average beam
        minor axis HPBW in deg), bpa (average beam position angle in
        deg) set the corresponding keywords in the headers. Choose
        '1420405745.51' for restfreq if you want the HI rest
        frequency.

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

            # Generate Gaussians or point sources and put properties in header

            for j in range(gauprops[i].shape[0]):
                if mode == 'gauss':
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
                else:
                    darray = np.zeros((naxis2,naxis1), dtype = target.dtype)
                    darray[int(gauprops[i][j,1]),int(gauprops[i][j,0])] = gauprops[i][j,2]
                    
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

    def gentestcubes(self, outcubi_prefix = ['test'], mode = 'gauss', naxis = 4, naxis1 = 257, naxis2 = 513, naxis3 = 4, pixelsize = 8.3333335E-04, ctype3 = 'VRAD', channelwidth = 5000., cellscal = '1/F', amp0 = 1., ainc = 0.0, pinc = 5., bmaj0 = 10., bmin0 = 7., binc = 0., bpa0 = 0, cinc = 0, avbeam = True):
        """Generate a set of test cubes with Gaussians or point sources

        Parameters:
        outcubi_prefix (list of str) : Prefixes of output cubes
        mode (str)                   : 'gauss' or 'point'
        naxis (int)                  : Number of axes 3 or 4
        naxis1 (int)                 : Dimension of cube axis 1
        naxis2 (int)                 : Dimension of cube axis 2
        naxis3 (int)                 : Dimension of cube axis 3 (# of channels)
        pixelsize (float)            : Pixelsize in degrees
        ctype3 (str)                 : Type of axis 3 'VRAD' or 'FREQ'
        channelwidth (float)         : Channel width (axis 3)
        cellscal (str)               : CELLSCAL keyword ('CONSTANT' or '1/F')
        amp0 (float)                 : Amplitude first Gaussian/point source
        ainc (float)                 : Increment in amplitude
        pinc (float)                 : Increment in position
        bmaj0 (float)                : HPBW major axis first Gaussian
        bmin0 (float)                : HPBW minor axis first Gaussian
        binc (float)                 : Increment in HPBWs
        bpa0 (float)                 : Beam position angle first Gaussian
        cinc (float)                 : Increment of beam position angles
        avbeam (bool)                : Provide average beam properties

        Creates n = len(outcubi_prefix) cubes with names
        outcubi_prefix[0]+'.fits', ..., outcubi_prefix[n-1]+'.fits'
        and sizes naxis1 x naxis2 x naxis3 and one additional fake
        'STOKES' axis 4. If mode == 'gauss', it will generate Gaussians, if
        'point' point sources (nearest pixel). The initial
        Gaussian/point source is centered at the centre of each plane,
        Then per plane the position changes along axis 1 and axis 2 by
        pinc pixels, and the amplitude by ainc. The first Gaussian is
        characterized by the HPBWs bmaj0, bmin0, and bpa0. bmaj and
        bmin incrase by binc per channel, bpa by cinc.

        """
        gauprops = []
        for i in range(len(outcubi_prefix)):
            planar = np.zeros((naxis3,6))
            for j in range(naxis3):
                planar[j, 0] = naxis1//2+(i*naxis3+j)*pinc
                planar[j, 1] = naxis2//2+(i*naxis3+j)*pinc
                planar[j, 2] = amp0+(i*naxis3+j)*ainc
                if mode == 'gauss':
                    planar[j, 3] = bmaj0+(i*naxis3+j)*binc
                    planar[j, 4] = bmin0+(i*naxis3+j)*binc
                    planar[j, 5] = bpa0+(i*naxis3+j)*cinc
            gauprops.append(planar)

        outcubi = [s + '.fits' for s in outcubi_prefix]
        if avbeam:
            bmaj = planar[:, 3].mean()*pixelsize
            bmin = planar[:, 4].mean()*pixelsize
            bpa = planar[:, 5].mean()
            
            self.createstcubes(mode = mode, gauprops = gauprops, outcubi = outcubi, naxis = naxis, naxis1 = naxis1, naxis2 = naxis2, ctype3 = ctype3, channelwidth = channelwidth, cellscal = cellscal, bmaj = bmaj, bmin = bmin, bpa = bpa)
        else:
            self.createstcubes(mode = mode, gauprops = gauprops, outcubi = outcubi, naxis = naxis, naxis1 = naxis1, naxis2 = naxis2, ctype3 = ctype3, channelwidth = channelwidth, cellscal = cellscal)


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

    for i in range(len(beach.binfo_input)):
        print('Input cube {}'.format(i))
        print('bmaj      : ', beach.binfo_input[i][:,0])
        print('bmin      : ', beach.binfo_input[i][:,1])
        print('bpa       : ', beach.binfo_input[i][:,2])
        print('freq      : ', beach.binfo_input[i][:,3])
        print('dx        : ', beach.binfo_input[i][:,4])
        print('rfreq/freq: ', beach.binfo_input[i][:,5])
        print('bmaj*freq : ', beach.binfo_input[i][:,6])
        print('bmin*freq : ', beach.binfo_input[i][:,7])
        print()
        print('Pixel cube {}'.format(i))
        print('bmaj      : ', beach.binfo_pixel[i][:,0])
        print('bmin      : ', beach.binfo_pixel[i][:,1])
        print('bpa       : ', beach.binfo_pixel[i][:,2])
        print('freq      : ', beach.binfo_pixel[i][:,3])
        print('dx        : ', beach.binfo_pixel[i][:,4])
        print('rfreq/freq: ', beach.binfo_pixel[i][:,5])
        print('bmaj*freq : ', beach.binfo_pixel[i][:,6])
        print('bmin*freq : ', beach.binfo_pixel[i][:,7])
        print()
        print('Target cube {}'.format(i))
        print('bmaj      : ', beach.binfo_target[i][:,0])
        print('bmin      : ', beach.binfo_target[i][:,1])
        print('bpa       : ', beach.binfo_target[i][:,2])
        print('freq      : ', beach.binfo_target[i][:,3])
        print('dx        : ', beach.binfo_target[i][:,4])
        print('rfreq/freq: ', beach.binfo_target[i][:,5])
        print('bmaj*freq : ', beach.binfo_target[i][:,6])
        print('bmin*freq : ', beach.binfo_target[i][:,7])
        print()

def gettoflux(incubi):
    """
    Derive total flux in a list of cubes
    """
    for incubus in incubi:
        hdul = fits.open(incubus)
        incubus_data = np.squeeze(hdul[0].data)
        incubus_header = hdul[0].header
        pixelsize = hdul[0].header['CDELT2']
        naxis3 = incubus_data.shape[0]

        output = {}

        for prefix in ['BMAJ', 'BMIN', 'BPA']:
            # Start with an empty array
            output[prefix] = np.empty((naxis3,))
            output[prefix][:] = np.nan
                    
            # Go through dict and search for a cubedefault
            cubedefault = np.nan
            for key in incubus_header.keys():
                if key == prefix:
                    cubedefault = incubus_header[prefix]
                    break

            # 'Channel' numbers
            chnum = []

            # Values
            chval = []
    
            # Go through dict and assign numbers
            for key in incubus_header.keys():
                if key.startswith(prefix):
                    if key != prefix:
                        try:
                            chnum = int(key[len(prefix):])
                            try:
                                chval = float(incubus_header[key])
                                output[prefix][chnum-1] = chval
                            except:
                                pass
                        except:
                            pass

            # Now fill remaining nans
            output[prefix][np.isnan(output[prefix])] = cubedefault
        for plane in range(naxis3):
            print('{}: Total flux plane {:d}: {:.3f}'.format(incubus, plane, incubus_data[plane,:].sum()*pixelsize*pixelsize/(2*np.pi/np.log(256.)*output['BMAJ'][plane]*output['BMIN'][plane])))

def test1():

    print()
    print('########################')
    print('########################')
    print(' Test 1: Frequency axis')
    print('########################')
    print()

    incubi = ['test1_incubus1', 'test1_incubus2']
    incubif = [i+'.fits' for i in incubi]
    Beach(verb = False).gentestcubes(outcubi_prefix = incubi, ctype3 = 'FREQ', channelwidth = 20000, amp0 = 1., ainc = 0.0, pinc = 5., bmaj0 = 10., bmin0 = 7., binc = 0., bpa0 = 0, cinc = 0)
    print('Test 1: Created input cubes {:s} and {:s}, with properties:'.format(incubif[0], incubif[1]))
    printcubeinfo(incubif[0])
    printcubeinfo(incubif[1])
    beach = Beach(inc_cubes = incubif, gentrans_exe = False)
    printbeachconts(beach)
    print()
    print('########################')

def test2():

    print()
    print('#################################')
    print('#################################')
    print(' Test 2: Velocity axis and model ')
    print('#################################')
    print()

    incubi =  ['test2_incubus1',  'test2_incubus2']
    incubif = [i+'.fits' for i in incubi]
    pincubi = ['test2_pincubus1', 'test2_pincubus2']
    pincubif = [i+'.fits' for i in pincubi]
    
    # 
    Beach(verb = False).gentestcubes(outcubi_prefix = incubi, amp0 = 1., ainc = 0.0, pinc = 5., bmaj0 = 10., bmin0 = 7., binc = 0., bpa0 = 0, cinc = 0, mode = 'gauss', ctype3='VRAD', channelwidth = 5000)
    
    print('Created input cubes {:s} and {:s}, with properties:'.format(incubi[0], incubi[1]))
    printcubeinfo(incubif[0])
    printcubeinfo(incubif[1])
    
    Beach(verb = False).gentestcubes(outcubi_prefix = pincubi, amp0 = 1., ainc = 0.0, pinc = -5., mode = 'point', ctype3='VRAD', channelwidth = 5000)

    beach = Beach(inc_cubes = incubif, gentrans_exe = False)

    outcubif = ['test2_outcubus1.fits', 'test2_outcubus2.fits']

    # This is just a preliminary exercise to create a parameter input interface
    beach = Beach(inc_cubes = incubif, gentrans_exe = False)
#    printbeachconts(beach)
    beach.gentrans(tra_fitsnames = outcubif, tra_modelnames = pincubif, tra_overwrite = True)
    print()
    print('Test 2: Created output cubes {:s} and {:s}'.format(outcubif[0], outcubif[1]))
    print()

    print()
    print('########################')

def test3():
    print('########################')
    print('########################')
    print('Test 3: Short version (shorter is still possible)')
    print('########################')
    print()

    incubi =  ['test3_incubus1',  'test3_incubus2']
    incubif = [i+'.fits' for i in incubi]
    pincubi = ['test3_pincubus1', 'test3_pincubus2']
    pincubif = [i+'.fits' for i in pincubi]
    
    Beach(verb = False).gentestcubes(outcubi_prefix = incubi, amp0 = 1., ainc = 0.0, pinc = 5., bmaj0 = 10., bmin0 = 7., binc = 0., bpa0 = 0, cinc = 0, mode = 'gauss', ctype3='VRAD', channelwidth = 5000)
    Beach(verb = False).gentestcubes(outcubi_prefix = pincubi, amp0 = 1., ainc = 0.0, pinc = -5., mode = 'point', ctype3='VRAD', channelwidth = 5000)

    # This is all it takes
    outcubif = ['test3_outcubus1.fits', 'test3_outcubus2.fits']
    Beach(inc_cubes = incubif, tra_modelnames = pincubif, tra_fitsnames = outcubif, tra_overwrite = True)
    print()
    print('Test 3: Created output cubes {:s} and {:s}'.format(outcubif[0], outcubif[1]))
    print()

def test4():

    print('################################')
    print('################################')
    print(' Test 4: Just images, not cubes')
    print('################################')
    print()
    
    incubi =  ['test4_incubus1',  'test4_incubus2']
    incubif = [i+'.fits' for i in incubi]
    pincubi = ['test4_pincubus1', 'test4_pincubus2']
    pincubif = [i+'.fits' for i in pincubi]
    
    Beach(verb = False).gentestcubes(outcubi_prefix = incubi, naxis = 2, naxis3 = 1, amp0 = 1., ainc = 0.0, pinc = 5., bmaj0 = 10., bmin0 = 7., binc = 0., bpa0 = 0, cinc = 0, mode = 'gauss', ctype3='VRAD', channelwidth = 5000)
    
    Beach(verb = False).gentestcubes(outcubi_prefix = pincubi, naxis = 2, naxis3 = 1, amp0 = 1., ainc = 0.0, pinc = -5., mode = 'point', ctype3='VRAD', channelwidth = 5000)

    print('Created input cubes {:s} and {:s}'.format(incubif[0], incubif[1]))
    #    printcubeinfo(incubi[0])
    #    printcubeinfo(incubi[1])
    outcubif = ['test4_outcubus1.fits', 'test4_outcubus2.fits']

    beach = Beach(inc_cubes = incubif, tra_modelnames = pincubif, tra_fitsnames = outcubif, gentrans_exe = False)
    printbeachconts(beach)
    print()
    print('Test 4: Created output cubes {:s} and {:s}'.format(outcubif[0], outcubif[1]))
    print()

def test5():
    print('########################')
    print('########################')
    print(' Test 5: Statistics plots')
    print('########################')
    print()
    incubi = []
    pincubi = []
    incubif = []
    pincubif = []
    outcubif = []
    for i in range(10):
        incubi.append('test5_incubus_{:d}'.format(i))
        pincubi.append('test5_pincubus_{:d}'.format(i))
        incubif.append('test5_incubus_{:d}.fits'.format(i))
        pincubif.append('test5_pincubus_{:d}.fits'.format(i))
        outcubif.append('test5_outcubus_{:d}.fits'.format(i))

    Beach(verb = False).gentestcubes(outcubi_prefix = incubi, naxis = 2, naxis3 = 20, amp0 = 1., ainc = 0.01, pinc = 0.5, bmaj0 = 10., bmin0 = 7., binc = 0.1, bpa0 = 0, cinc = 0.3, mode = 'gauss', ctype3='VRAD', channelwidth = 5000)
    
    Beach(verb = False).gentestcubes(outcubi_prefix = pincubi, naxis = 2, naxis3 = 20, amp0 = 1., ainc = 0.01, pinc = -0.5, bmaj0 = 10., cinc = 0.3, mode = 'point', ctype3='VRAD', channelwidth = 5000)
    
    beach = Beach(inc_cubes = incubif, gentrans_exe = False, bst_percents = 95)
    printbeachconts(beach)
    beach.genhistoplots(hist_plotname='test5.png', hist_interactive='test5.html', hist_scaling = 'constant', hist_sample = 'chan', hist_n_per_bin = 2, hist_overwrite = True)
    #beach.gentrans(tra_fitsnames = outcubif, tra_modelnames = pincubif, tra_overwrite = True)
    print()
    print('Test 5: Created output plots')
    print()
    
def test6():
    print('########################')
    print('########################')
    print(' Test 6: Masking')
    print('########################')
    print()
    incubi = []
    pincubi = []
    incubif = []
    pincubif = []
    outcubif_mask = []
    outcubif_scale = []
    outcubif_hybrid = []
    outcubif_max = []
    for i in range(1):
        incubi.append('test6_incubus_{:d}'.format(i))
        pincubi.append('test6_pincubus_{:d}'.format(i))
        incubif.append('test6_incubus_{:d}.fits'.format(i))
        pincubif.append('test6_pincubus_{:d}.fits'.format(i))
        outcubif_scale.append('test6_outcubus_scale_{:d}.fits'.format(i))
        outcubif_mask.append('test6_outcubus_mask_{:d}.fits'.format(i))
        outcubif_hybrid.append('test6_outcubus_hybrid_{:d}.fits'.format(i))
        outcubif_max.append('test6_outcubus_max_{:d}.fits'.format(i))

    Beach(verb = False).gentestcubes(outcubi_prefix = incubi, naxis = 2, naxis3 = 20, amp0 = 1., ainc = 0.00, pinc = 0.5, bmaj0 = 10., bmin0 = 7., binc = 0.5, bpa0 = 0, cinc = 0.3, mode = 'gauss', ctype3='VRAD', channelwidth = 5000)
    
    #Beach(verb = False).gentestcubes(outcubi_prefix = pincubi, naxis = 2, naxis3 = 20, amp0 = 1., ainc = 0.01, pinc = -0.5, bmaj0 = 10., cinc = 0.3, mode = 'point', ctype3='VRAD', channelwidth = 5000)
    
    beach = Beach(inc_cubes = incubif, gentrans_exe = False, bst_percents = 50)
    printbeachconts(beach)
    #beach.genhistoplots(hist_plotname='test6.png', hist_interactive='test6.html', hist_scaling = 'constant', hist_sample = 'total', hist_n_per_bin = 2, hist_overwrite = True)
    #beach.gentrans(tra_fitsnames = outcubif, tra_modelnames = pincubif, tra_overwrite = True)
    beach.gentrans(tra_fitsnames = outcubif_mask, tra_mode = 'mask', tra_tol = 0.001, tra_maxker = 1E7, tra_overwrite = True)
    beach.gentrans(tra_fitsnames = outcubif_scale, tra_mode = 'scale', tra_tol = 0.001, tra_maxker = 1E7, tra_overwrite = True)
    beach.gentrans(tra_fitsnames = outcubif_hybrid, tra_mode = 'hybrid', tra_tol = 0.001, tra_maxker = 1E7, tra_overwrite = True)
    beach.gentrans(tra_fitsnames = outcubif_max, tra_mode = 'max', tra_tol = 0.001, tra_maxker = 1E7, tra_overwrite = True)
    print()
    print('Test 6: Created output cubes')
    print()
    print('Reading output cubes and determining approximate total flux per channel')
    gettoflux(incubif)
    gettoflux(outcubif_mask)
    gettoflux(outcubif_scale)
    gettoflux(outcubif_hybrid)
    gettoflux(outcubif_max)
    
def testing():
    
    #Beach().convoltests()

    print('')
    print('##################')
    print('##################')
    print('##################')
    print(' Final Tests/Demo ')
    print('##################')
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    
def description():
    """
    Describing the class's properties
    """
    return '\n' + \
textwrap.fill('This module is a radioastronomical tool. Its purpose is to transform a set of images with known resolution, which may vary from image to image, into a set of images with the same resolution (by choice relative to a frequency reference frame, see below):')+'\n'+\
'\n'+\
textwrap.fill('The resolution of a radioastronomical image is usually represented by a two-dimensional Gaussian, the (clean) beam, whose properties are known, and which is described by having an amplitude of 1 and a major- and minor axis of the ellipse at its half-power level (major axis \"half-power-beam-width\" (HPBW), or minor axis HPBW), as well as the position angle of the major axis (measuered anticlockwise from the North). The intensity (spectral brightness) in radioastronomical images is assumed to be the true sky brightness convolved with the individual Gaussians. In the following, parameters or quantities with suffixes, mediumfixes, or prefixes \'bmaj\' or \'BMAJ\' are related to the beam major axis HPBW, parameters or quantities with suffixes, mediumfixes, or prefixes \'bmin\' or \'BMIN\' are related to the beam major axis HPBW, parameters or quantities with suffixes, mediumfixes, or prefixes \'bpa\' or \'BPA\' are related to the beam majore axis position angle.')+'\n'+\
'\n'+\
textwrap.fill('The module takes two sets of images or spectroscopic data cubes in FITS format as an input. One is assumed to contain (part of) the sky brightness convolved with a beam (a \"restored image\" or a \"resiudual\"), the other is assumed to be a sky model, which is not yet convolved with the beam. Equolver re-convolves each image/plane in the first data set to a common beam, or a set of common beams, which can be shared among all cubes, all planes, or within individual cubes. whose properties can be derived from the statistics of the known beams. To do so, the images are Fourier-transformed, then divided by the Fourier-transform of the original beam, multiplied with the Fourier-transform of the target beam, and Fourier-transformed back. Alternatively the images are scaled with the integral of the target beam (the target beam-solid-angle, BSA) divided by the BSA of the original beam. Hybrid approaches are also possible, see below. The images/planes in the second data set (the model) only get convolved with the target beam. Then both images are added.')+'\n'+\
'\n'+\
textwrap.fill('In the following we provide a detailed description of the module and its input parameters. In some cases, the manual refers to lists as input parameters. Those are entered in Python style, as comma-separated lists of values enclosed by opening and closing brackets. Strings not entered inside a list can be entered without quotes, inside lists, strings should be enclosed by single- and double quotes. As many shell interpreters would interpret quotes themselves, the user has to take care for the right format. If the user wants e.g. to enter the list [\'foo\', \'fooly\'], on the command line in bash or csh this may have to be entered as e.g. --parameter \"[\'foo\', \'fooly\']\" (notice the double quotes).')+'\n'+\
'\n'+\
textwrap.fill('Equolver potentially goes through four steps, with the parameter prefixes indicating to which step they belong:')+'\n'+\
'\n'+\
textwrap.fill('- Reading cubes (gen\'inc\'ubus)', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('- Reading the beam information (gen\'bin\'fo)', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('- Generating beam statistics (gen\'bst\'ats)', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('- Generating (interactive) histograms for diagnostics (gen\'hist\'oplots)', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('- Generating (a) common beam(s) from statistics or direct input (gen\'tar\'get)', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('- Generate a set of transformed images from the input with the beam properties derived (gen\'tra\'ns)', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
'\n'+\
'Reading cubes (gen\'inc\'ubus):'+'\n'+\
textwrap.fill('============================')+'\n'+\
'\n'+\
textwrap.fill('The user specifies the input cubes meant to be re-convolved with --inc_cubes INC_CUBES or -i INC_CUBES, where INC_CUBES is either a string with the name of the single data cube (image) or a list of strings with the names of the cubes (images).')+'\n'+\
'\n'+\
textwrap.fill('Expansion scheme for parameters:')+'\n'+\
textwrap.fill('===============================')+'\n'+\
'\n'+\
textwrap.fill('Some parameters can take multiple formats, referring to the input list of cubes/images: The user can enter single values for all planes in all cubes/images. The user can instead enter lists, where each value of the list corresponds to all planes in a data cube with the same index. If the user enters a list of lists, each member of the list (which is also a list) corresponds to a data cube/image with the same index, and its elements correspond to the planes in the corresponding data cube. Values can be entered in astropy syntax as strings of as floats, where the units are degrees or Hz.')+'\n'+\
'\n'+\
textwrap.fill('Example:')+'\n'+\
'\n'+\
textwrap.fill('cube1.fits and cube2.fits each contain 3 planes (with indices from 0 to 2), and the user specifies -i \"[\'cube1\', \'cube2\']\".')+'\n'+\
'\n'+\
textwrap.fill('- The user specifies --bin_bmaj 0.002 -> for all planes in both cubes the input bmaj is 0.002 degrees', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('- The user specifies --bin_bmaj \"[0.002, \'1 arcmin\']\" -> for all planes in cube1.fits the input bmaj is 0.002 degrees, for all planes in cube2.fits the input bmaj is 1 arcminute.', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('- The user specifies --bin_bmaj \"[[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]]\" -> input bmaj is 0.01 deg for cube1.fits plane 0, 0.02 deg for cube1.fits plane 1, 0.03 deg for cube1.fits plane 2, 0.04 deg for cube2.fits plane 0, 0.05 deg for cube2.fits plane 1, 0.06 deg for cube2.fits plane 2.', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
'\n'+\
textwrap.fill('Reading the beam information (gen\'bin\'fo):')+'\n'+\
textwrap.fill('=========================================')+'\n'+\
'\n'+\
textwrap.fill('If present, the beam properties of the individual images (planes) are read from the FITS headers of INC_CUBES (see \'Reading cubes\'), in which they have a format BMAJ_NNN, BMIN_NNN, BPA_NNN, where NNN is the plane number, starting with 1. If the keywords BMAJ, BMIN, BPA (without suffix) are present in the header, they serve as a default value for the respective plane-specific specifications. Alternatively, the user can directly provide the default values using --bin_bmaj BIN_BMAJ, --bin_bmin BIN_BMIN, and --bin_bpa BIN_BPA, following the expansion scheme for parameters: using this direct input, the user can provide one number for all channels and cubes/images, a list of numbers, providing one number per cube, and a list of lists, providing a list of numbers (per channel) for each cube. If --bin_bmaj_replace True is set, the header values BMAJ or BMAJ_NNN in the data cubes are ignored and the input (\"default\") values are read in instead. If --bin_bmin_replace True is set, the header values BMIN or BMIN_NNN in the data cubes are ignored and the input (\"default\") values are read in instead. If --bin_bpa_replace True is set, the header values BPA or BPA_NNN in the data cubes are ignored and the input (\"default\") values are read in instead.')+'\n'+\
'\n'+\
textwrap.fill('By nature, the third axis of a radiointerferometric data cube is either frequency, or velocity. Equolver interprets any velocity as radio velocity VRAD with respect to a rest frequency nu0: VRAD = c*(nu0-nu)/nu0, where c is the speed of light and nu the frequency in a specific channel. The rest frequency nu0 is read from the header of each cube using the keyword \'RESTFREQ\'. Defaults can be entered using the parameter --bin_restfreq BIN_RESTFREQ, which can be entered using the expansion scheme, and enforcing the usage of this direct input can be achieved by setting --bin_restfreq_replace True .')+'\n'+\
'\n'+\
textwrap.fill('In a radio observation, without adding any additional weighting scheme for the visibilities that changes with frequency, the beam major and minor axes scale with the inverse of the frequency. With equolver, the user has the possibility to scale the beam to a \'normalisation\' frequency before deriving statistics and/or scaling it back to the specific frequency per channel. This way, a common beam for the normalisation frequency can be calculated to then scale this beam to the specific channels, to then reconvolve the cube planes. With the parameter --bin_normfreq BIN_NORMFREQ (defaulting to 1.4 GHz) the user can enter this normalisation frequency. With the normalisation-frequency nf (see above), assuming that the beam sizes scale with 1/frequency we derive b(nf) = b(f)*f/nf , where b is major or minor axis HPBW. Also this parameter is expandable, but we strongly recommend not to use more than one value for all data cubes.')+'\n'+\
'\n'+\
textwrap.fill('Generating beam statistics (gen\'bst\'ats):')+'\n'+\
textwrap.fill('========================================')+'\n'+\
'\n'+\
textwrap.fill('After reading the beam information, the user can generate statistics on the collected beam properties. The section can be omitted by setting --genbstats_exe False (Notice that this disables consecutive sections except plotting). The parameters --bst_parameter BST_PARAMETER, --bst_scaling BST_SCALING, --bst_sample BST_SAMPLE, and --bst_stype BST_STYPE determine which statistics are being calculated (by default all statistics are calculated, but a choice can accelerate the processing time). If for any of the parameters \'all\' is chosen (which is the default), all fields are filled. If the scaling type is \'constant\', the given (read) values for bmaj and bmin are evaluated, if \'frequency\' is chosen, all beam sizes are scaled to the same normalisation-frequency nf (see above), assuming that the beam sizes scale with 1/frequency. b(nf) = b(f)*f/nf . In the following, bsa is the beam solid angle of an individual (or average) beam, ceb the circular equivalent beam, the circular beam with a beam solid angle of bsa (sqrt(bmaj*bmin)). The parameters can be combined, and the following can be calculated:')+'\n'+\
'\n'+\
textwrap.fill('--bst_parameter BST_PARAMETER', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('BST_PARAMETER:', initial_indent=8*' ', subsequent_indent=10*' ')+'\n'+\
textwrap.fill('\'bmaj\'       major axis dispersions/hpbws', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'bmin\'       minor axis dispersions/hpbws', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'bpa\'        beam position angles',         initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'bsa\'        beam solid angle',             initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'ceb\'        circular equivalent beam dispersions/hpbws', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
'\n'+\
textwrap.fill('--bst_scaling BST_SCALING', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('BST_SCALING:', initial_indent=8*' ', subsequent_indent=10*' ')+'\n'+\
textwrap.fill('\'const\'      constant', initial_indent='            ', subsequent_indent='              ')+'\n'+\
textwrap.fill('\'frequency\'  frequency', initial_indent='            ', subsequent_indent='              ')+'\n'+\
'\n'+\
textwrap.fill('--bst_stype BST_STYPE', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('BST_STYPE:', initial_indent=8*' ', subsequent_indent=10*' ')+'\n'+\
textwrap.fill('\'minimum\'    Minimum',                   initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'maximum\'    Maximum',                   initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'average\'    Average',                   initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'stdev\'      Standard deviation',        initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'median\'     Median',                    initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'mad\'        Median-absolute-deviation', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'madstdev\'   Standard deviation calculated from the median-absolute-deviation', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'percentile\' Score at percents',         initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'commonbeam\' Common beam as calculated using the radio-beam module (https://radio-beam.readthedocs.io) based on the Khachiyan algorithm (https://en.wikipedia.org/wiki/Ellipsoid_method). Parameters bst_tolerance, bst_nsamps, epsilon are used for this method', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
'\n'+\
textwrap.fill('--bst_sample BST_SAMPLE', initial_indent=4*' ', subsequent_indent=6*' ')+'\n'+\
textwrap.fill('BST_SAMPLE:', initial_indent=8*' ', subsequent_indent=10*' ')+'\n'+\
textwrap.fill('\'cube\'       Statistics to be carried out for all channels per cube (generates lists with length of the number of input cubes)', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'chan\'       Statistics to be carried out for all cubes per channel (generates lists with a length of the maximum number of channels in any cube)', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
textwrap.fill('\'total\'      Statistics to be carried out for all channels in all cubes (generates a float)', initial_indent=12*' ', subsequent_indent=25*' ')+'\n'+\
'\n'+\
textwrap.fill('For some of those parameters, additional arguments are required. --bst_percents BST_PERCENTS is the number of percents for the percentile statistics, while --bst_tolerance BST_TOLERANCE, the convergence tolerance, --bst_nsamps BST_NSAMPS, number of vertices of the polygon describing the elliptic beam, and --bst_epsilon, the allowance to overestimate the common beam area, --bst_maxiter BST_MAXITER, the maximum number of iterations to find a common beam, are specific to the calculation of a common beam (see above, see https://radio-beam.readthedocs.io). Notice that equolver will not stop to search for a common beam. If a search fails, it will be re-started with a higher number of iterations and a higher tolerance.')+'\n'+\
'\n'+\
textwrap.fill('Generating (interactive) histograms for diagnostics (gen\'hist\'oplots)')+'\n'+\
textwrap.fill('=====================================================================')+'\n'+\
'\n'+\
textwrap.fill('Generates histograms of the beam properties and the statistics read in and generated in the previous sections. This will either produce static histograms (in png format) with the name HIST_PLOTNAME if the parameter --hist_plotname HIST_PLOTNAME is set and/or an interactive html file with the name HIST_INTERACTIVE if the parameter --hist_interactive HIST_INTERACTIVE is set. With the parameter --hist_sample HIST_SAMPLE the user can choose which statistics should be shown: \'cube\' means that the histogram is generated for each cube, \'chan\' means that it is generated for a specific channel in all cubes, \'total\' (the default) means that the statistics for all channels in all cubes are shown in one plot. The parameter --hist_scaling HIST_SCALING decides if the statistics are plotted using the original beam properties (if HIST_SCALING is set to \'constant\') or if they are first scaled to the norm frequency (if HIST_SCALING is set to \'frequency\') assuming that the beam scales proportionally to the inverse of the frequency. With --hist_overwrite HIST_OVERWRITE the user indicates whether existing files can be overwritten (if HIST_OVERWRITE is set to \'True\', the default is False).')+'\n'+\
'\n'+\
textwrap.fill('Generating (a) common beam(s) from statistics or direct input (gen\'tar\'get)')+'\n'+\
textwrap.fill(70*'=')+'\n'+\
'\n'+\
textwrap.fill('The section generates the target beam properties (bmaj, bmin, bpa) for all input data cubes and all planes therein. The section can be omitted by setting --gentarget_exe False (Notice that this disables the consecutive section). They are generated as follows: For each quantity bmaj, bmin, bpa the parameters')+'\n'+\
'\n'+\
textwrap.fill('--tar_quant_inter TAR_QUANT_INTER (intercept, default: 0)', initial_indent=4*' ', subsequent_indent=22*' ')+'\n'+\
textwrap.fill('--tar_quant_slope TAR_QUANT_SLOPE (slope, default: 0)', initial_indent=4*' ', subsequent_indent=22*' ')+'\n'+\
textwrap.fill('--tar_quant_absc  TAR_QUANT_ABSC (abscissa, default: 0)', initial_indent=4*' ', subsequent_indent=22*' ')+'\n'+\
'\n'+\
textwrap.fill('result in the output quantity quant calculated as:')+'\n'+\
'\n'+\
textwrap.fill('quant = TAR_QUANT_INTER + TAR_QUANT_SLOPE*TAR_QUANT_ABSC', initial_indent=4*' ', subsequent_indent=13*' ')+'\n'+\
'\n'+\
textwrap.fill('The parameters are either direct inputs of the target quantities following the expansion scheme described above. Alternatively, TAR_QUANT_INTER and TAR_QUANT_ABSC can be a list of four strings, denoting (in that order) bst_parameter, bst_scaling, bst_stype, and bst_sample as described in the gen\'bst\'ats section above. The corresponding values from the statistics will then be used. Notice that \'freq\' means the inverse of the transformation of the input beam info struct into the frequency-scaled one. Major and minor axis beams entered directly are then interpreted as beams at a frequency as specified in the parameter --normfreq NORMFREQ.')+'\n'+\
'\n'+\
textwrap.fill('Example (slightly artificial):')+'\n'+\
'\n'+\
textwrap.fill('--tar_bmaj_inter \"[\'bmaj\', \'frequency\', \'percentile\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_slope 0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_absc  0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_inter \"[\'bmaj\', \'frequency\', \'average\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_slope \"[\'bmaj\', \'frequency\', \'stdev\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_absc  3.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_inter  \"[\'bpa\', \'frequency\',  \'average\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_slope  0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_absc   0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_scaling    frequency'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
'\n'+\
textwrap.fill('This would use the percentile (at BST_PERCENTS) of the major axis (scaled to the normalisation frequency) of all planes in all data cubes for the target major axis, the average of the major axis HPBWs (scaled to the normalisation frequency) of all planes in all data cubes plus 3 times the standard deviation of all planes in all data cubes for the target minor axes, and the average position angle of all planes in all data cubes. All results are scaled back from the normalisation frequency to the actual frequency in each channel.')+'\n'+\
'\n'+\
'Sensible choices are:''\n'+\
'\n'+\
textwrap.fill('Use the (90%) percentile of the major axis HPWBs at the normalisation frequency as the target major and minor axes, scale back with frequency. This is the default:')+'\n'+\
'\n'+\
textwrap.fill('--tar_bmaj_inter \"[\'bmaj\', \'frequency\', \'percentile\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_slope 0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_absc  0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_inter \"[\'bmaj\', \'frequency\',  \'percentile\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_slope 0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_absc  0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_inter  0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_slope  0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_absc   0.0'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_scaling    frequency'+'\n', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
'\n'+\
textwrap.fill('Use the average of the major axis HPWBs  at the normalisation frequency plus 3 sigma for major and minor HPBW, scale with 1/F:')+'\n'+\
'\n'+\
textwrap.fill('--tar_bmaj_inter \"[\'bmaj\', \'frequency\', \'average\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_slope \"[\'bmaj\', \'frequency\', \'stdev\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_absc  3.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_inter \"[\'bmaj\', \'frequency\',  \'percentile\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_slope \"[\'bmaj\', \'frequency\', \'stdev\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_absc  3.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
'\n'+\
textwrap.fill('Use the common beam at the normalisation frequency, scale with 1/F:')+'\n'+\
'\n'+\
textwrap.fill('--tar_bmaj_inter \"[\'bmaj\', \'frequency\', \'commonbeam\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_slope 0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmaj_absc  0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_inter \"[\'bmin\', \'frequency\',  \'commonbeam\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_slope 0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bmin_absc  0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_inter  \"[\'bpa\', \'frequency\',  \'commonbeam\', \'total\']\"', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_slope  0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
textwrap.fill('--tar_bpa_absc   0.0', initial_indent=4*' ', subsequent_indent=21*' ')+'\n'+\
'\n'+\
textwrap.fill('Generate a set of transformed images from the input with the beam properties derived (gen\'tra\'ns)')+'\n'+\
textwrap.fill(70*'=')+'\n'+\
'\n'+\
textwrap.fill('The final section serially opens all cubes (images) listed in inc_cubes and generates cubes (de-)convolved to the resolution as specified in the target section. The section can be omitted by setting --gentrans_exe False. If executed the section optionally convolves all cubes listed in --tra_modelnames TRA_MODELNAMES with the respective Gaussians and adds those to the output. If given, the number of cubes listed as TRA_MODELNAMES and the dimensionality of the cubes must be identical to the one of the input cubes (specified through INC_CUBES). While from a mathematical viewpoint the de-convolution to smaller beams should work, this is in practice limited by numerical effects. To handle such situations, equolver provides several strategies, specified with the parameter --tra_mode TRA_MODE. A (de-)convolution is declared a success or a failure by comparing for the input and output plane the sum of the pixels of the plane divided by the beam solid angle. If the ratio of the larger sum divided by the smaller sum is larger than the parameter --tra_tol TRA_TOL (default: 2), the total power hence changes significally, the (de-)convolution is flagged a failure. In addition, --tra_maxker TRA_MAXKER allows the user to manually set the maximum in the FT of a kernel as an alternative way to decide whether a reconvolution is bound to fail.')+'\n'+\
'\n'+\
'The mode of the deconvolution is determined with TRA_MODE as follows:'+'\n'+\
'\n'+\
textwrap.fill('\'scale\':  Do not convolve but scale the intensity to the target beam (divide by original beam solid angle and multiply with target beam solid angle.)', initial_indent=4*' ', subsequent_indent=14*' ')+'\n'+\
textwrap.fill('\'mask\':   (De-)convolve and mask channel in the output if the deconvolution fails.', initial_indent=4*' ', subsequent_indent=14*' ')+'\n'+\
textwrap.fill('\'hybrid\': Attempt to (de-) convolve the plane and fall back to \'scale\' if the (de-)convolution fails.', initial_indent=4*' ', subsequent_indent=14*' ')+'\n'+\
textwrap.fill('\'max\':    Attempt to (de-) convolve the plane. If that fails, convolve along the beam minor axis to the target minor beam if that is larger than the original. Then scale.', initial_indent=4*' ', subsequent_indent=14*' ')+'\n'+\
'\n'+\
textwrap.fill('If parameter --tra_hdmode TRA_HDMODE is set to True, a keyword \'EQMODE\' with the value of TRA_MODE is added to the header. In addition, if TRA_MODE is not set to \'scale\' or \'mask\', any plane for which the first convolution failed, is highlighted by the keyword-value pair EQS_i = \'SCALE\' or EQSC = \'HYBRID\'. i in this context is the plane number (Fortran/FITS style: starting with 1).')+'\n'+\
'\n'+\
textwrap.fill('Finally, the user can generate beam properties in the header of the output cubes. If --tra_commonbeam TRA_COMMONBEAM is set to \'True\', the average beam properties are calculated and inserted into the header as \'BMAJ\', \'BMIN\', \'BPA\'. In addition the keyword \'BEAMSCAL\' is added to the header and its value set to \'CONSTANT\' if the beam frequency scaling is de-activated and hence TAR_SCALING is set to \'False\'. It is set to \'1/F\' if the beam frequency scaling is activated and hence TAR_SCALING is set to \'True\'.')+\
'\n'

def parsing():
    if '-v' in sys.argv or '--verb' in sys.argv:
        epilog = description()
    else:
        epilog = 'Use \'equolver -h -v\' for verbose description.'
    parser = argparse.ArgumentParser(description='Convolve fits images and data cubes to the same resolution.', formatter_class=argparse.RawTextHelpFormatter, prog='equolver', usage='%(prog)s [options]', epilog = epilog, fromfile_prefix_chars='@', argument_default=argparse.SUPPRESS)

    # Common
    parser.add_argument('--threads', help='Number of threads')
    parser.add_argument('--version', action = 'version', version = version)
    parser.add_argument('--verb', '-v',    help='Verbose output? \'True\': yes')
    
    # gen
    parser.add_argument('--inc_cubes', '-i', help='Input cubes: names or list of names, python style')

    # genbinput
    parser.add_argument('--bin_bmaj',             help='Beam major axis default value(s), format see below')
    parser.add_argument('--bin_bmaj_replace',     help='Enforce usage of default values bin_bmaj? (True = yes)')
    parser.add_argument('--bin_bmin',             help='Beam minor axis default value(s), format see below')
    parser.add_argument('--bin_bmin_replace',     help='Enforce usage of default values bin_bmin? (True = yes)')
    parser.add_argument('--bin_bpa',              help='Beam position angle default value(s), format see below')
    parser.add_argument('--bin_bpa_replace',      help='Enforce usage of default values bin_bpa? (True = yes)')
    parser.add_argument('--bin_restfreq',         help='Rest frequency default value(s), format see below')
    parser.add_argument('--bin_restfreq_replace', help='Enforce usage of default values bin_restfreq?')
    parser.add_argument('--bin_normfreq',         help='Frequency in Hz to normalize beam to if mode is \'frequency\'')
    
    parser.add_argument('--genbstats_exe', help='Generate beam statistics (\'True\': yes)')
    parser.add_argument('--bst_parameter', help='Parameter name (\'all\', \'bmaj\', \'bmin\', \'bpa\', \'bsa\', \'ceb\')')
    parser.add_argument('--bst_scaling',   help='Scaling type (\'all\', \'constant\', \'frequency\')')
    parser.add_argument('--bst_stype',     help='Type of statistics to calculate (\'all\', \'minimum\', \'maximum\', \'average\', \'stdev\', \'median\', \'mad\', \'madstdev\', \'percentile\', \'percents\', \'commonbeam\')')
    parser.add_argument('--bst_sample',    help='Sample(s) to calculate statistics on (\'all\', \'cube\', \'chan\', \'total\')')
    parser.add_argument('--bst_percents',  help='Percents for the percentile statistics')
    parser.add_argument('--bst_tolerance', help='Tolerance for searching the common beam')
    parser.add_argument('--bst_nsamps',    help='Number of edges of beam for common beam')
    parser.add_argument('--bst_epsilon',   help='Epsilon to search for common beam')
    parser.add_argument('--bst_maxiter',   help='Maximum iterations to search for common beam')

    parser.add_argument('--hist_plotname',    help='Name of static plot')
    parser.add_argument('--hist_interactive', help='Name of interactive plot')
    parser.add_argument('--hist_sample',      help='Sample to plot \'cube\', \'chan\', or \'total\'')
    parser.add_argument('--hist_scaling',     help='Scaling to use (\'frequency\' or \'constant\')')
    parser.add_argument('--hist_overwrite',   help='Allow overwriting files produced before')

    parser.add_argument('--gentarget_exe',     help='Generate target beams (\'True\': yes)')
    parser.add_argument('--tar_bmaj_inter',    help='Beam major axis intercept')
    parser.add_argument('--tar_bmaj_slope',    help='Beam major axis slope')
    parser.add_argument('--tar_bmaj_absc',     help='Beam major axis abscissae')
    parser.add_argument('--tar_bmin_inter',    help='Beam minor axis intercept')
    parser.add_argument('--tar_bmin_slope',    help='Beam minor axis slope')
    parser.add_argument('--tar_bmin_absc',     help='Beam minor axis abscissae')
    parser.add_argument('--tar_bpa_inter',     help='Beam position angle axis intercept')
    parser.add_argument('--tar_bpa_slope',     help='Beam position angle axis slope')
    parser.add_argument('--tar_bpa_absc',      help='Beam position angle abscissae')
    parser.add_argument('--tar_scaling',       help='Use 1/F scaling when calculating the target array, either \'frequency\' or \'input\'')

    parser.add_argument('--gentrans_exe',   help='Generate transformed cubes (\'True\': yes)')
    parser.add_argument('--tra_modelnames', '-m', help='Input fits file names, containing the models.')
    parser.add_argument('--tra_fitsnames', '-o', help='Output fits file names.')
    parser.add_argument('--tra_mode',       help='\'scale\', \'mask\', \'hybrid\', \'max\'')
    parser.add_argument('--tra_tol',        help='tolerance to determine if convolution failed')
    parser.add_argument('--tra_maxker',     help='Maximum value that the FT of the convolution kernel can assume, will assume failure if larger')
    parser.add_argument('--tra_commonbeam', help='Generate common (average) beam information in header')
    parser.add_argument('--tra_indibeam',   help='Generate individual beam information in header')
    parser.add_argument('--tra_hdmode',     help='Generate information about scaling/convolution in header')
    parser.add_argument('--tra_overwrite',  help='Overwrite output if already existent (\'True\': yes)?')

    whatnot = parser.parse_args()
    inpars = vars(whatnot)

    for key in inpars.keys():
        try:
            result = eval(inpars[key])
        except:
            result = inpars[key]
        inpars[key] = result

    # if 'inc_cubes' in inpars.keys():
    #     print(inpars['inc_cubes'])
    #     if inpars['inc_cubes'] == True:
    #         print('yo')
    
    return inpars
    
def runtime():
    kwargs = parsing()
    for argument in ['help', 'version']:
        if argument in kwargs.keys():
            sys.exit()
    Beach(**kwargs)    

if __name__ == '__main__':
    #testing()
    runtime()
