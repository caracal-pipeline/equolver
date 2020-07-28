#!/usr/bin/env python
import warnings
import copy
import numpy as np
from astropy.io import fits
from astropy import constants
from astropy import wcs
from astropy import units
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
                       bmaj = None, bmaj_replace = None,
                       bmin = None, bmin_replace = None,
                       bpa = None, bpa_replace = None,
                       restfreq = HIFREQ, restfreq_replace = None,
                       loadcubes = True):
        """
        Private instance variables:
        (multiple: None, a float, a list of floats, a numpy array, 
         or a list with numpy arrays)

        _cubenames (list)     : Names of input cubes
        _cubes (list)         : List of open cube hdus

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
                                tion, in the order of cubenames and
                                cubes slowest index same as index of
                                cube in cubes, followed by _bmaj,
                                _bmin, and _bpa, fastest index
                                channels, so it is an array of size
                                (len(_cubenames),3,chans), where
                                chans is the number of channels.
        _binfo_pixel (list)  : _binfo_input converted into pixel
                               scaling using dispersion instead of HPBW
        _binfo_pixel_f (list): _beampixelprops divided by frequency
                               (where appropriate) in GHz
        _bstats (dict)       : Dictionary containing all statistics
        _bgoal_pixel_f (list): Target properties of beams per plane
                               per plane per cube
        _bgoal_input (list)  : Target beam properties
        """
        self.resetvars()
        self.initcubes(cubenames = cubenames)
        self.initbinfo_input()
        
        self.getdefault(bmaj, 'bmaj')
        self.getdefault(bmin, 'bmin')
        self.getdefault(bpa, 'bpa')
        self.getdefault(restfreq, 'restfreq')

        self.setdefreplace(bmaj_replace = bmaj_replace,
                           bmin_replace = bmin_replace,
                           bpa_replace = bpa_replace,
                           restfreq_replace= restfreq_replace)

        if loadcubes:
            if self._cubenames != None:
                self.loadcubes()
        return

    def resetvars(self):
        """
        Reset/init instance variables
        """

        # Cubes
        self._cubenames = None
        self._cubes = None

        # Defaults
        self._bmaj = None
        self._bmaj_replace = None
        self._bmin = None
        self._bmin_replace = None
        self._bpa = None
        self._bpa_replace = None
        self._restfreq = None
        self._restfreq_replace = None

        # bmaj, bmin, pa, nu, pixelsize
        self._binfo_input = None

        # bmaj, bmin, sin pa, cos pa, nu
        self._binfo_pixel = None
        
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
        self.initcubes(cubenames = value)
        return
    
    @cubenames.deleter
    def cubenames(self, value):
        self.resetcubes()
        
    @property
    def cubes(self):
        """
        Return a copy of the list cubes

        The cubes themselves are not copies
        """
        # Notice that this is not a deep copy to save memory
        # But this is really the responsibility of the user
        
        if type(self._cubes) == type(None):
            return None
        return copy.copy(self._cubes)

    def returndefault(self, value):
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
        return self.returndefault(self._bmaj)

    @bmaj.setter
    def bmaj(self, value):
        """
        Set bmaj
        """
        self.getdefault(bmaj, 'bmaj')
        if value != self._bmaj:
            warning('Setting bmaj has failed. Likely, you need to load cubes',
                    'first')
        return

    @bmaj.deleter
    def bmaj(self, value):
        self._bmaj = None
        return
        
    @property
    def bmin(self):
        """
        Return a copy of bmin
        """
        return self.returndefault(self._bmin)

    @bmin.setter
    def bmin(self, value):
        """
        Set bmin
        """
        self.getdefault(bmin, 'bmin')
        if value != self._bmin:
            warning('Setting bmin has failed. Likely, you need to load cubes',
                    'first')
        return

    @bmin.deleter
    def bmin(self, value):
        self._bmin = None
        return

    @property
    def bpa(self):
        """
        Return a copy of rchan
        """
        return self.returndefault(self._bpa)

    @bpa.setter
    def bpa(self, value):
        """
        Set bpa
        """
        self.getdefault(bpa, 'bpa')
        if value != self._bpa:
            warning('Setting bpa has failed. Likely, you need to load cubes',
                    'first')
        return

    @bpa.deleter
    def bpa(self, value):
        self._bpa = None
        return

    @property
    def restfreq(self):
        """
        Return a copy of restfreq
        """
        return self.returndefault(self._restfreq)

    @restfreq.setter
    def restfreq(self, value):
        """
        Set restfreq
        """
        self.getdefault(restfreq, 'restfreq')
        if value != self._restfreq:
            warning('Setting restfreq has failed. Likely, you need to load cubes first')
        return

    @restfreq.deleter
    def restfreq(self, value):
        self._restfreq = None
        return

    @property
    def bmaj_replace(self):
        return self._bmaj_replace

    @bmaj_replace.setter
    def bmaj_replace(self, value):
        self._bmaj_replace = value
    
    @bmaj_replace.deleter
    def bmaj_replace(self, value):
        self._bmaj_replace = None
    
    @property
    def bmin_replace(self):
        return self._bmaj_replace

    @bmin_replace.setter
    def bmin_replace(self, value):
        self._bmin_replace = value
    
    @bmin_replace.deleter
    def bmin_replace(self, value):
        self._bmin_replace = None
    
    @property
    def bpa_replace(self):
        return self._bpa_replace

    @bpa_replace.setter
    def bpa_replace(self, value):
        self._bpa_replace = value
    
    @bpa_replace.deleter
    def bpa_replace(self, value):
        self._bpa_replace = None
    
    @property
    def restfreq_replace(self):
        return self._restfreq_replace

    @restfreq_replace.setter
    def restfreq_replace(self, value):
        self._restfreq_replace = value
        return
    
    @restfreq_replace.deleter
    def restfreq_replace(self, value):
        self._bmaj_replace = None
    
    @property
    def binfo_input(self):
        """
        Return a copy of binfo_input
        """
        return self.returndefault(self._binfo_input)

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
        return self.returndefault(self._binfo_pixel)

    @property
    def binfo_pixel_f(self):
        """
        Return a copy of binfo_pixel_
        """
        return self.returndefault(self._binfo_pixel_f)

    @property
    def bstats(self):
        """
        Return a copy of bstats
        """
        return copy.deepcopy(self._bstats)

    def resetcubes(self):
        """
        Close all cubes in instance and set cubes and cubenames to None
        """
        if self._cubes != None:
            for cube in self._cubes:
                cube.close()
            self._cubes = None
        self._cubenames = None
        return

    def initcubes(self, cubenames = None):
        """
        Add cubenames to intrinsic cubenames

        Input:
        cubenames (list of str): List of names of input cubes/images

        Reads cubenames in as the list of target cubes/images
        """
        if type(cubenames) == type(None):
            return
        self.resetcubes()
        self._cubenames = cubenames[:]
        self._cubes = []
        for cube in self._cubenames:
            self._cubes += [fits.open(cube)]
        return
    
    def initbinfo_input(self):
        """
        Initiate the beam info arrays for the cubes and fill with nans

        Will initiate a list of ndarrays, one for each cube in 
        self.cubes with the dimension nchan x 5, where nchan is the
        number of channels in the corresponding cube. Each value is
        set to np.nan
        """
        
        if type(self._cubes) == type(None):
            self.initcubes()
            
        if type(self._cubes) == type(None):
            return

        self._binfo_input = []
        
        for cube in self._cubes:
            
            # Create array of nans, size naxis3 x 5 (bmaj, bmin, bpa,
            # psize, freq) and append to the info list
            if 'NAXIS3' in cube[0].header.keys():
                naxis3 = cube[0].header['NAXIS3']
            else:
                naxis3 = 1
        

            self._binfo_input.append(np.empty((naxis3,5,)))
            self._binfo_input[-1][:] = np.nan
        return

    def getdefault(self, inquant, quantname = None):
        """
        Fill target quantname, wich is the name of an instance
        variable, with the content of inquant
        
        Input:
        (multiple: None, a float, a list of floats, a numpy array, 
         or a list with numpy arrays)
        inquant (multiple): input
        quantname (str)   : Name of quantity (without leading _)

        self._quantname is assumed to be a list of linear ndarrays of
        the size of naxis3 of self._cubes (i.e. the number of
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
        to have as many elements as self._cubes and an ndarray always
        has to have the same dimension as the naxis3 of the
        corresponding cube.

        """
        if type(inquant) == type(None) and \
           type(self.__dict__['_'+quantname]) == type(None):
            inquant = np.nan
        
        # Check if there is something sensible asked for
        if type(quantname) == type(None):
            return
        
        if type(inquant) == type(None):
            return

        # Invoke initbinfo_input
        if type(self._binfo_input) == type(None):
            self.initbinfo_input()

        # If binfo_input is empty, return
        if type(self._binfo_input) == type(None):
            return

        self.__dict__['_'+quantname] = []

        for i in range(len(self._cubes)):

            # One number for all
            if type(inquant) == type(.1) or type(inquant) == type(np.nan):
                self.__dict__['_'+quantname].append(np.repeat(
                    inquant, self._binfo_input[i].shape[0]))

            # If input is an ndarray it is a default array with
            # one inquant per channel for all cubes
            elif type(inquant) == type(np.array([])): 
                self.__dict__['_'+quantname].append(np.repeat(
                    inquant, self._binfo_input[i].shape[0]))

            # If it is a list, there are two possibilities
            elif type(inquant) == type([]):

                # One number for all channels
                if type(inquant[i]) == type(.1) or \
                type(inquant[i]) == type(np.nan):
                    self.__dict__['_'+quantname].append(np.repeat(
                        inquant[i], self._binfo_input[i].shape[0]))

                # The only other possibility is an ndarray
                elif type(inquant[i]) == type(np.array([])):
                    self.__dict__['_'+quantname].append(inquant[i])

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

        return

    def setdefreplace(self, bmaj_replace = None,
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
            
    def genbinfo_input(self, bmaj = None, bmaj_replace = None,
                          bmin = None, bmin_replace = None,
                          bpa = None,  bpa_replace = None,
                          restfreq = HIFREQ, restfreq_replace = None):
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

        if type(self._cubes) == type(None):
            warnings.warn('No cubes loaded.')
            return
        
        if type(self._binfo_input) == type(None):
            self.initbinfo_input()
            
        if type(self._binfo_input) == type(None):
            warning('binfo_input not present. Use initcubes first.')
            return
        
        # Get defaults
        self.getdefault(bmaj, 'bmaj')
        self.getdefault(bmin, 'bmin')
        self.getdefault(bpa, 'bpa')
        self.getdefault(restfreq, 'restfreq')

        self.setdefreplace(bmaj_replace = bmaj_replace,
                           bmin_replace = bmin_replace,
                           bpa_replace = bpa_replace,
                           restfreq_replace = restfreq_replace)

        for i in range(len(self._cubes)):

            # It would be rather surprising if there is a channel-dependent
            # rest frequency, but who knows
            restfreqtab = self.getchanval(
                'RESTFREQ', self._cubes[i][0].header, value = self._restfreq[i],
                usevalue = restfreq_replace, usedefault = True)
            
            # Count numbers of axes
            naxis = self._cubes[i][0].header['NAXIS']

            pixcrd = np.column_stack((
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64),
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64)))
            pixcrd2 = np.column_stack((
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64),
                np.ones(self._binfo_input[i].shape[0], dtype=np.float64)+1.))

            stocol = 3
            if naxis > 2:
                if self._cubes[i][0].header['CTYPE3'] != 'STOKES':
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
                    self._cubes[i][0].header['RESTFREQ'] = restfreqtab[j]
                wcshand = kwcs.Projection(self._cubes[i][0].header)
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
            if 'CELLSCAL' in self._cubes[i][0].header.keys():
                if self._cubes[i][0].header['CELLSCAL'] == '1/F':
                    if not type(wcshand2) == type(None):
                        if not np.isnan(self._binfo_input[i][:,3].sum()):
                            self._binfo_input[i][:,4] = wcshand2.crval[2]* \
                                self._binfo_input[i][:,4]/ \
                                self._binfo_input[i][:,3]

            cellscal_use_constant = True
            if type(wcshand2) == type(None):
                cellscal_use_constant = False
            if 'CELLSCAL' in self._cubes[i][0].header.keys():
                if self._cubes[i][0].header['CELLSCAL'] == '1/F':
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
            
            # finlist: slowest index cube in 'cubes' list, followed
            # by bmaj, bmin, and bpa
            self._binfo_input[i][:,0] = self.getchanval( \
                'BMAJ', self._cubes[i][0].header, value = self._bmaj[i],
                usevalue = bmaj_replace, usedefault = True, dscal = dscal)
            self._binfo_input[i][:,1] = self.getchanval( \
                'BMIN', self._cubes[i][0].header, value = self._bmin[i],
                usevalue = bmin_replace, usedefault = True, dscal = dscal)
            self._binfo_input[i][:,2] = self.getchanval( \
                'BPA', self._cubes[i][0].header, value = self._bpa[i],
                usevalue = bpa_replace, usedefault = True)
        return

    def getchanval(self, prefix, thedict, value = None, usevalue = False,
                   usedefault = True, dscal = 1.):
        """Return an ndarray of values from an input dict with formatted
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
        getchanval('BMAJ', ourdict, value = 17.3, usevalue = False, 
                   usedefault = True)
        [0.1, 0.3, 0.2]

        ourdict = { 'BMAJ1': 0.1, 'BMAJ3': 0.2, 'NAXIS3': 3}
        getchanval('BMAJ', ourdict, value = 17.3, usevalue = False,
                   usedefault = True)
        [0.1, 17.3, 0.2]

        getchanval('BMAJ', ourdict, value = [17.2,17.3,17.4], 
                   usevalue = True, usedefault = True)
        [17.2, 17.3, 17.4]

        getchanval('BMAJ', ourdict, value = None, usevalue = False, 
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

    def genbinfo_pixel(self):
        """
        Generate binfo_pixel, a list of ndarrays, which is
        the beam info file binfo_input converted to pixel
        units.
        """

        
        if type(self._binfo_input) == type(None):
            warnings.warn('Trying to generate beam info.')
            self.genbinfo_input()
            
        if type(self._binfo_input) == type(None):
            warnings.warn('No beam information read, which disables',
                          ' further processing.')
            return
        
        for binfarray in self._binfo_input:
            if binfarray[:,0].sum() == np.nan:
                warnings.warn('Not sufficient information about beam major ',
                              'axes, which disables further processing.')
                return
        
            if binfarray[:,1].sum() == np.nan:
                warnings.warn('Not sufficient information about beam minor ',
                              'axes, which disables further processing.')
                return
            
            if binfarray[:,2].sum() == np.nan:
                warnings.warn('Not sufficient information about beam position ',
                              'angle, which disables further processing.')
                return
            
            if binfarray[:,4].sum() == np.nan:
                warnings.warn('Not sufficient information about pixel ',
                              'dimension, which disables further processing.')
                return
            
        self._binfo_pixel = []
        for binfarray in self._binfo_input:
            boutfarray = binfarray.copy()
            
            # Scale to pixels as unit and dispersion
            boutfarray[:,0] = binfarray[:,0]/binfarray[:,4]/np.sqrt(np.log(256))
            boutfarray[:,1] = binfarray[:,1]/binfarray[:,4]/np.sqrt(np.log(256))

            # Scale to cos pa and sin pa
            #boutfarray[:,2] = np.sin(np.pi*binfarray[:,2]/180.)
            #boutfarray[:,3] = np.cos(np.pi*binfarray[:,2]/180.)

            # Copy frequency
            # boutfarray[:,4] = binfarray[:,3]
            
            self._binfo_pixel.append(boutfarray)
        exit
        
    def genbinfo_pixel_f(self):
        """
        Generate binfo_pixel, a list of ndarrays, which is
        the beam info file binfo_input converted to pixel
        units, divided by frequency.
        """
        if type(self._binfo_pixel) == type(None):
            warnings.warn('Trying to generate beam info.')
            self.genbinfo_pixel()

        if type(self._binfo_pixel) == type(None):
            warnings.warn('No beam information read in, which disables',
                          ' further processing.')
            return
        
        self._binfo_pixel_f = []
        for bixfarray in self._binfo_pixel:
            bouxfarray = bixfarray.copy()
            if not np.isnan(bixfarray[:,3].sum()):
                print('got here2')
                bouxfarray[:,0] = bixfarray[:,0]/(bixfarray[:,3]/1E9)
                bouxfarray[:,1] = bixfarray[:,1]/(bixfarray[:,3]/1E9)
            self._binfo_pixel_f.append(bouxfarray)

    def initbstats(self):
        """Create a structure able to contain all statistics

        The intrinsic statistics structure _bstats contains all
        statistics. It is a nested dictionary. See genbstats for a
        description. This method generates the structure and fills it
        with nans.

        """

        # Count the number of cubes and the maximum number of channels
        if type(self._binfo_input) == type(None):
            warnings.warn('Trying to generate beam info.')
            self.genbinfo_input()

        if type(self._binfo_input) == type(None):
            warnings.warn('No beam information read in, which disables',
                          ' further processing.')
            return

        ncubes = len(self._binfo_input)

        channels = 0
        for i in range(ncubes):
            channels = max(self._binfo_input[i].shape[0], channels)

        if channels == 0:
            warnings.warn('No cube with more than 0 channels present.'+\
                          ' This is virtually impossible. No statistics '+\
                          'structure built.')
            return
            
        self._bstats = {}
        for key0 in ['bmaj', 'bmin', 'pang']:
            self._bstats[key0] = {}
            for key1 in ['input', 'pixel', 'pixel_f']:
                self._bstats[key0][key1] = {}
                for key2 in ['minimum', 'maximum', 'average', 'stddev', 'median',
                             'mad', 'madstddev', 'percentile', 'percentage',
                             'commonbeam']:
                    self._bstats[key0][key1][key2] = {}
                    # for key4 in ['cube', 'chan', 'tota']:
                    self._bstats[key0][key1][key2]['total'] = np.nan
                    self._bstats[key0][key1][key2]['chan'] = \
                        np.empty((channels,))
                    self._bstats[key0][key1][key2]['chan'][:] = np.nan
                    self._bstats[key0][key1][key2]['cube'] = []
                    for i in range(ncubes):
                        self._bstats[key0][key1][key2]['cube'].append(np.nan)
                        
        return
                    
    def genbstats(self, parameter = 'all', scaling = 'all',
                  stype = 'all', perc = 90, sample = 'all',):
        """
        Generate statistics and dump it into the bstats structure
        
        Input:
        parameter (str or list of str): Parameter name ('bmaj', 'bmin', 'pang')
        scaling (str or list of str)  : Scaling type ('input', 'pixel',
                                        'pixel_f')
        stype (str or list of str)    : Type of statistics to
                                        calculate ('minimum',
                                        'maximum', 'average',
                                        'stdev', 'median', 'mad',
                                        'medstdev', 'percentile',
                                        'percent', 'combeam')
        perc (float)                  : Percentile for the percentile
                                        statistics
        sample (str or list of str)   : Sample(s) to calculate
                                        statistics on ('cube', 'chan',
                                        'total')

        The method generates statistics on the collected beam
        properties. The parameters parameter, scaling, sample, and
        type determine which part of the bstats structure gets
        filled. If for any of the parameters 'all' is chosen (which is
        the default), all fields are filled. The scaling type
        corresponds to the tables/structures binfo_input, binfo_pixel,
        and binfo_pixel_f generated using methods genbinfo_input,
        genbinfo_pixel, genbinfo_pixel_f.

        From top to bottom level:

        parameter: 
            major axis dispersions/hpbws ('bmaj')
            minor axis dispersions/hpbws ('bmin')
            beam position angles ('pang')

        scaling (corresponding to generated structs):
            Input (beam HPBWs in deg, PA in deg) ('input')
            Pixel (beam dispersions in pixel units, PA in deg) ('pixel')
            Pixel per frequency (beam dispersions in pixel units,
            divided by Frequency in GHz, PA in deg) ('pixel_f')

        stype:
            Minimum ('minimum')
            Maximum ('maximum')
            Average ('average')
            Standard deviation ('stddev')
            Median ('median')
            Median-absolute-deviation ('mad')
            Standard deviation calculated from the
                median-absolute-deviation ('madstddev')
            x%-percentile ('percentile')
            Percentage percentile is corresponding to ('percentage')

            Common beam as calculated using the radio-beam module
            (https://radio-beam.readthedocs.io) based on the Khachiyan
            algorithm (https://en.wikipedia.org/wiki/Ellipsoid_method)
                ('commonbeam')

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

        # First expand all parameters to the same format
        if type(parameter) == type(''):
            parameter = [parameter]
        if type(scaling) == type(''):
            scaling = [scaling]
        if type(stype) == type(''):
            stype = [stype]
        if type(sample) == type(''):
            sample = [sample]

        para = parameter[:]
        scal = scaling[:]
        styp = stype[:]
        samp = sample[:]

        if 'all' in para:
            para = ['bmaj', 'bmin', 'pang']
        if 'all' in scal:
            scal = ['input', 'pixel', 'pixel_f']
        if 'all' in styp:
            styp = ['minimum', 'maximum', 'average', 'stddev', 'median',
                    'mad', 'madstddev', 'percentile', 'percentage',
                    'commonbeam']
        if 'all' in samp:
            samp = ['cube', 'chan', 'total']
        
        # Make sure we do have the required information available
        for item in scale:
            if type(self.__dict__['_binfo_'+quantname]) == type(None):
                method = getattr(self, 'genbinfo_'+quantname)
                warnings.warn('Attempting to genereate {:s}'.format(item),
                              ' information struct.')
                method()
            if type(self.__dict__['_binfo_'+quantname]) == type(None):
                warnings.warn('{:s} information struct cannot be generated. ',
                              'Returning without generating statistics.')
                return

        # Now let's do it
        for sca in scal:

            # get struct
            struct = self.__dict__['_binfo_'+sca])
            
            for sam in samp:
                
                # collect sample as list of np arrays
                collist = []
                if sam == 'cube':
                    for i in struct:
                        collist.append[i[:,[0,1,2]]]
                elif sam == 'chan':
                    for i in range(self._bstats[sca][sam][bmaj]['chan']):
                        # Continue here
                        pass
                elif sam == 'total':
                    pass
                else:
                    warnings.warn('Statistics cannot be generated. Check sample '
                                  'parameter. Returning without generating ',
                                  'statistics.')
                    
                for sty in styp:
                    # Get function to generate statistics
                    for par in para:
                        # Apply statistics
                        # Put the right value back
                        pass
                    
    _genminimum(self, parray, )
    
    def loadcubes(self, cubenames = None,
                  bmaj = None, bmaj_replace = None,
                  bmin = None, bmin_replace = None,
                  bpa = None, bpa_replace = None,
                  restfreq = None, restfreq_replace = None):
        """
        Load a list of cubes and relevant quantities from their 
        headers
        
        Input:
        cubenames (list of str): List of input cube names
        bmaj (None or list):     Beam major axis default value(s)
        bmaj_replace (bool):     Enforce usage of default values?
                                 (True = yes)
        bmin (None or list):     Beam minor axis default value(s)
        bmin_replace (bool):     Enforce usage of default values?
                                 (True = yes)
        bpa (None or list):      Beam position angle default value(s)
        bpa_replace (bool):      Enforce usage of default values? 
                                 (True = yes)
        restfreq (None or list): Rest frequency default value(s)
        restfreq_replace (bool): Enforce usage of default values?
                                 (True = yes, defaults to False)        
        """
        if type(restfreq) == type(None):
            restfrequency = self.HIFREQ

        # Read in cubes
        self.initcubes(cubenames = cubenames)

        # Get the beam properties for every plane in every cube
        self.genbinfo_input(bmaj = bmaj, bmin = bmin, bpa = bpa,
                             restfreq = restfreq, bmaj_replace = bmaj_replace,
                             bmin_replace = bmin_replace,
                             bpa_replace = bpa_replace,
                             restfreq_replace= restfreq_replace)
        self.genbinfo_pixel()
        self.genbinfo_pixel_f()
        self.initbstats()
        
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
    print('Pixel cube 1 freq')
    print('bmaj    : ', beach.binfo_pixel_f[0][:,0])
    print('bmin    : ', beach.binfo_pixel_f[0][:,1])
    print('cos bpa : ', beach.binfo_pixel_f[0][:,2])
    print('sin bpa : ', beach.binfo_pixel_f[0][:,3])
    print('freq    : ', beach.binfo_pixel_f[0][:,4])
    print('Pixel cube 2 freq')
    print('bmaj    : ', beach.binfo_pixel_f[1][:,0])
    print('bmin    : ', beach.binfo_pixel_f[1][:,1])
    print('cos bpa : ', beach.binfo_pixel_f[1][:,2])
    print('sin bpa : ', beach.binfo_pixel_f[1][:,3])
    print('freq    : ', beach.binfo_pixel_f[1][:,4])
    print()
    print(beach.bstats)

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
    beach.loadcubes(cubenames = params['cubes'],restfreq=beach.HIFREQ)
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
    beach = Beach()
    beach.loadcubes(cubenames = params['cubes'],restfreq=beach.HIFREQ)
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
    beach = Beach()
    beach.loadcubes(cubenames = params['cubes'],restfreq=beach.HIFREQ)
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
    beach = Beach()
    beach.loadcubes(cubenames = params['cubes'],restfreq=beach.HIFREQ)
    printbeachconts(beach)
