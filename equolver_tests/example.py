#!/usr/bin/env python
import numpy as np
from astropy.io import fits

def getchanval(prefix, thedict, default = 0.0):
    """Return an ndarray of values from an input dict with formatted keywords

    Input:
    prefix (string): Prefix of formatted keywords
    thedict (dict) : Input dict
    default (float): Default value
    
    Rather specialized method. Assume that there is a keyword 'pref'
    number of keywords 'prefix1', 'prefix2', ..., in the list. Then
    the method returns an array of values sorted from prefix1 to
    prefixn where n is the maximum number after prefix. Any number not
    found gets a value assigned to prefix, if that does not exist as a
    keyword, default.

    Example: 
    ourdict = { 'BMAJ1': 0.1, 'BMAJ3': 0.2, 'BMAJ': 0.3}
    getchanval('BMAJ', ourdict, default = 17.3)
        [0.1 0.3 0.2]

    ourdict = { 'BMAJ1': 0.1, 'BMAJ3': 0.2}
    getchanval('BMAJ', ourdict, default = 17.3)
        [0.1 17.3 0.2]
    """
    
    # 'Channel' numbers
    chnum = np.empty((0))

    # Values
    chval = np.empty((0))

    # Go through dict and search for prefix, then append to both arrays
    searchdefault = True
    for key in thedict.keys():
        if key.startswith(prefix):
            if key[len(prefix):] == '':
                if searchdefault:
                    thedefault = thedict[key]
                    searchdefault = False
            else:
                chnum = np.append(chnum, np.array([int(key[len(prefix):])]).reshape(1), axis=0)
                chval = np.append(chval, np.array([thedict[key]]).reshape(1), axis=0)
                
    # Get max channel and fill the output array, using defaults if required
    finlist = np.empty((0,1))
    for i in range(int(chnum.max())):
        hulp = np.nonzero(chnum == i+1)[0]
        if hulp.size > 0:
            index = int(hulp[0])
            finlist = np.append(finlist, chval[int(hulp[0])])
        else:
            finlist = np.append(finlist, thedefault)

    # Finis
    return finlist
    
params = { 'cubes': ['mh_p23_HI_t10_r0_1370_1390.image.fits',
                     'mh_p24_HI_t10_r0_1370_1390.image.fits']
}

for cube in params['cubes']:
    hdu = fits.open(cube)
    finlist = getchanval('BMAJ', hdu[0].header)
    hdu.close()
    
