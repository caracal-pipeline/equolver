#!/usr/bin/env python
import warnings
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
                       restfreq = HIFREQ, restfreq_replace = None):
        """
        Private instance variables:
        (multiple: None, a float, a list of floats, a numpy array, 
         or a list with numpy arrays)

        _cubenames (list)       : Names of input cubes
        _cubes (list)           : List of open cube hdus

        _bmaj     (multiple)    : Beam major axis default value(s)
        _bmaj_replace (bool)    : Enforce usage of default values?
                                  (True = yes)
        _bmin     (multiple)    : Beam minor axis default value(s)
        _bmin_replace (bool)    : Enforce usage of default values?
                                  (True = yes)
        _bpa     (multiple)     : Beam position angle default value(s)
        _bpa_replace (bool)     : Enforce usage of default values? 

        _restfreq (multiple)    : Rest frequency default value(s)
        _bpa_replace (bool)     : Enforce usage of default values? 

        _beamproplist (ndarray) : List of arrays of point spread func-
                                  tion, in the order of cubenames and
                                  cubes slowest index same as index of
                                  cube in cubes, followed by _bmaj,
                                  _bmin, and _bpa, fastest index
                                  channels, so it is an array of size
                                  (len(_cubenames),3,chans), where
                                  chans is the number of channels.
        """
        self.resetvars()
        self.initcubes(cubenames = cubenames)
        self.initbeamproplist()

        # Nonstandard header items are allowed a default
        self.getdefault(bmaj, 'bmaj')
        self.getdefault(bmin, 'bmin')
        self.getdefault(bpa, 'bpa')
        self.getdefault(restfreq, 'restfreq')
        
        self._bmaj_replace = False
        self._bmin_replace = False
        self._bpa_replace = False
        self._restfreq_replace = False

        if bmaj_replace != None:
            self._bmaj_replace = bmaj_replace
        if bmin_replace != None:
            self._bmin_replace = bmin_replace
        if bpa_replace != None:
            self._bpa_replace = bpa_replace
        if restfreq_replace != None:
            self._restfreq_replace = restfreq_replace
            
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

        # Central input array
        self._beamproplist = None

        return

    @property
    def cubenames(self):
        """
        Return a copy of cubenames
        """
        if type(self._cubenames) == type(None):
            return None
        return copy.deepcopy(self._beamproplist)

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
    def beamproplist(self):
        """
        Return a copy of beamproplist
        """
        return self.returndefault(self._beamproplist)

    @beamproplist.setter
    def beamproplist(self, value):
        """
        Copy the input
        """
        if type(value) == type(None):
            self._beamproplist = None
            return

        # Input must be a list of np arrays, which we will copy
        returnlist = []
        for i in value:
            returnlist.append(i.copy())
        self._beamproplist = returnlist
        return

    @beamproplist.deleter
    def beamproplist(self, value):
        self._beamproplist = None
    
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
    
    def initbeamproplist(self):
        """
        Generate the beam info arrays for the cubes and fill with nans

        Will generate a list of ndarrays, one for each cube in 
        self.cubes with the dimension nchan x 5, where nchan is the
        number of channels in the corresponding cube. Each value is
        set to np.nan
        """
        
        if type(self._cubes) == type(None):
            self.initcubes()
            
        if type(self._cubes) == type(None):
            return

        self._beamproplist = []
        
        for cube in self._cubes:
            
            # Create array of nans, size naxis3 x 5 (bmaj, bmin, bpa,
            # psize, freq) and append to the info list
            if 'NAXIS3' in cube[0].header.keys():
                naxis3 = cube[0].header['NAXIS3']
            else:
                naxis3 = 1
        

            self._beamproplist.append(np.empty((naxis3,5,)))
            self._beamproplist[-1][:] = np.nan
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

        if type(quantname) == type(None) and \
           type(self.__dict__['_'+quantname]) == type(None):
            quantname = np.nan
        
        # Check if there is something sensible asked for
        if type(quantname) == type(None):
            return
        
        if type(inquant) == type(None):
            return

        # Invoke initbeamproplist
        if type(self._beamproplist) == type(None):
            self.initbeamproplist()

        # If beamproplist is empty, return
        if type(self._beamproplist) == type(None):
            return

        self.__dict__['_'+quantname] = []

        for i in range(len(self._cubes)):

            # One number for all
            if type(inquant) == type(.1) or type(inquant) == type(np.nan):
                self.__dict__['_'+quantname].append(np.repeat(
                    inquant, self._beamproplist[i].shape[0]))

            # If input is an ndarray it is a default array with
            # one inquant per channel for all cubes
            elif type(inquant) == type(np.array([])): 
                self.__dict__['_'+quantname].append(np.repeat(
                    inquant, self._beamproplist[i].shape[0]))

            # If it is a list, there are two possibilities
            elif type(inquant) == type([]):

                # One number for all channels
                if type(inquant[i]) == type(.1) or \
                type(inquant[i]) == type(np.nan):
                    self.__dict__['_'+quantname].append(np.repeat(
                        inquant[i], self._beamproplist[i].shape[0]))

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

    def getbeamprop(self, bmaj = None, bmaj_replace = None,
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
            warning('No cubes loaded. Use initcubes first.')
            return
        
        if type(self._beamproplist) == type(None):
            self.initbeamproplist()
            
        if type(self._beamproplist) == type(None):
            warning('beamproplist not present. Use initcubes first.')
            return
        
        # Get defaults
        self.getdefault(bmaj, 'bmaj')
        self.getdefault(bmin, 'bmin')
        self.getdefault(bpa, 'bpa')
        self.getdefault(restfreq, 'restfreq')

        # We could do the same with replace but that's too much
        if type(self._bmaj_replace) == type(None) or \
                                   type(bmaj_replace) != type(None):
            self._bmaj_replace = bmaj_replace
        if type(self._bmin_replace) == type(None) or \
                                   type(bmin_replace) != type(None):
            self._bmin_replace = bmin_replace
        if type(self._bpa_replace) == type(None) or \
                                   type(bpa_replace) != type(None):
            self._bpa_replace = bpa_replace
        if type(self._restfreq_replace) == type(None) or \
                                   type(restfreq_replace) != type(None):
            self._restfreq_replace = restfreq_replace

        for i in range(len(self._cubes)):
            # finlist: slowest index cube in 'cubes' list, followed
            # by bmaj, bmin, and bpa
            self._beamproplist[i][:,0] = self.getchanval( \
                'BMAJ', self._cubes[i][0].header, value = self._bmaj[i],
                usevalue = bmaj_replace, usedefault = True)
            self._beamproplist[i][:,1] = self.getchanval( \
                'BMIN', self._cubes[i][0].header, value = self._bmin[i],
                usevalue = bmin_replace, usedefault = True)
            self._beamproplist[i][:,2] = self.getchanval( \
                'BPA', self._cubes[i][0].header, value = self._bpa[i],
                usevalue = bpa_replace, usedefault = True)

            # It would be rather surprising if there is a channel-dependent
            # rest frequency, but who knows
            restfreqtab = self.getchanval(
                'RESTFREQ', self._cubes[i][0].header, value = self._restfreq[i],
                usevalue = restfreq_replace, usedefault = True)
            
            # Count numbers of axes
            naxis = self._cubes[i][0].header['NAXIS']

            pixcrd = np.column_stack((
                np.ones(self._beamproplist[i].shape[0], dtype=np.float64),
                np.ones(self._beamproplist[i].shape[0], dtype=np.float64)))
            pixcrd2 = np.column_stack((
                np.ones(self._beamproplist[i].shape[0], dtype=np.float64),
                np.ones(self._beamproplist[i].shape[0], dtype=np.float64)+1.))

            stocol = 3
            if naxis > 2:
                if self._cubes[i][0].header['CTYPE3'] != 'STOKES':
                    pixcrd = np.column_stack((
                        pixcrd, np.arange(self._beamproplist[i].shape[0],
                                          dtype=np.float64)+1.))
                    pixcrd2 = np.column_stack((
                        pixcrd2, np.arange(self._beamproplist[i].shape[0],
                                           dtype=np.float64)+1.))
                else:
                    stocol = 2
                
            for j in range(stocol,naxis):             
                pixcrd = np.column_stack((
                    pixcrd, np.ones(self._beamproplist[i].shape[0],
                                    dtype=np.float64)))
                pixcrd2 = np.column_stack((
                    pixcrd2, np.ones(self._beamproplist[i].shape[0],
                                     dtype=np.float64)))

            # Use Kapteyn to find out about frequency and pixel size
            worldcrd = []
            worldcrd2 = []
            wcshand2 = None
            for j in range(self._beamproplist[i].shape[0]):
                if not np.isnan(restfreqtab[j]):
                    self._cubes[i][0].header['RESTFREQ'] = restfreqtab[j]
                wcshand = kwcs.Projection(self._cubes[i][0].header)
                if stocol == 3:
                    wcshand2 = wcshand.spectra('FREQ')
                    worldcrd += [wcshand2.toworld(pixcrd[j])[2]]
                worldcrd2 += [np.fabs(
                    wcshand2.toworld(pixcrd[j])[1]-
                    wcshand2.toworld(pixcrd2[j])[1])]

            if worldcrd != []:
                self._beamproplist[i][:,3] = np.array(worldcrd)
            self._beamproplist[i][:,4] = np.array(worldcrd2)

            # Now check for cellscal and change the pixel size
            # accordingly

            if 'CELLSCAL' in self._cubes[i][0].header.keys():
                if not type(wcshand2) == type(None):
                    if not np.isnan(self._beamproplist[i][:,3].sum()):
                        self._beamproplist[i][:,4] = wcshand2['CRVAL3']* \
                                   self._beamproplist[i][:,4]/ \
                                   self._beamproplist[i][:,3]
        return

    def getchanval(self, prefix, thedict, value = None, usevalue = False,
                   usedefault = True):
        """
        Return an ndarray of values from an input dict with formatted
        keywords

        Input:
        prefix (string)         : Prefix of formatted keywords
        thedict (dict)          : Input dict
        value (ndarray)         : Default values
        usevalue (None or bool) : If set output is identical to
                                  value (or np.nan if value is None)
        
    
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
        order not to use it as a default).

        Example: 
        ourdict = { 'BMAJ1': 0.1, 'BMAJ3': 0.2, 'BMAJ': 0.3, 'NAXIS3': 3}
        getchanval('BMAJ', ourdict, value = 17.3, usevalue = False, usedefault = True)
        [0.1, 0.3, 0.2]

        ourdict = { 'BMAJ1': 0.1, 'BMAJ3': 0.2, 'NAXIS3': 3}
        getchanval('BMAJ', ourdict, value = 17.3, usevalue = False, usedefault = True)
        [0.1, 17.3, 0.2]

        getchanval('BMAJ', ourdict, value = [17.2,17.3,17.4], usevalue = True, usedefault = True)
        [17.2, 17.3, 17.4]

        getchanval('BMAJ', ourdict, value = None, usevalue = False, usedefault = False)
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

        if type(value) != type(None):
            if type(value) == type(np.array([])):
                output[np.isnan(output)] = value[np.isnan(output)]
            else:
                output[np.isnan(output)] = value
        
        return output

    def topix(self):
        """
        Take array of beam properties in cube and translate to pixel coordinates and rad (beam)
        """
        pass
        
    def loadcubes(self, cubenames = None,
                  bmaj = None, bmaj_replace = False,
                  bmin = None, bmin_replace = False,
                  bpa = None, bpa_replace = False,
                  psize = None, psize_replace = False,
                  rchan = None, rchan_replace = False,
                  rfreq = None, rfreq_replace = False,
                  dfreq = None, dfreq_replace = False,
                  restfrequency = None,
                  restfrequency_replace = False):
        """
        Load a list of cubes and relevant quantities from their headers
        
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
        psize (None or list):    Pixel size default value(s)
        psize_replace (bool):    Enforce usage of default values?
                                 (True = yes)
        rchan (None or list):    Reference channel default value(s)
        rchan_replace (bool):    Enforce usage of default values?
                                 (True = yes)
        rfreq (None or list):    Reference frequency default value(s)
        rfreq_replace (bool):    Enforce usage of default values?
                                 (True = yes)
        dfreq (None or list):    Channel separation default value(s)
        dfreq_replace (bool):    Enforce usage of default values?
                                 (True = yes)
        restfreq (None or list): Rest frequency default value(s)
        restfreq_replace (bool): Enforce usage of default values?
                                 (True = yes)        
        """
        if type(restfrequency) == type(None):
            restfrequency = self.HIFREQ

        # Read in cubes
        self.initcubes(cubenames = cubenames)

        # Generate the beamproplist arrays
        self.initbeamproplist()

        # Get the beam properties for every plane in every cube
        self.getbeamprop(bmaj = bmaj, bmaj_replace = bmaj_replace,
                         bmin = bmin, bmin_replace = bmin_replace,
                         bpa = bpa, bpa_replace = bpa_replace)
        
        return
        
        


if __name__ == '__main__':
    params = { 'cubes': ['p1_02.fits',
                     'p2_02.fits']
    }
    beach = Beach()
    beach.loadcubes(cubenames = params['cubes'])
    print(beach.beamproplist[0][:,4])
    params = { 'cubes': ['p1_04.fits',
                     'p2_04.fits']
    }
    beach = Beach()
    beach.loadcubes(cubenames = params['cubes'],restfrequency=beach.HIFREQ)
    print(beach.beamproplist[1][:,4])
    

