#!/usr/bin/env python
import numpy as np

from equolver import beach
print('')
print('##############')
print(' Beach Demo ')
print('##############')

# The following parameters are used to create two data cubes using
# method createstcubes

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
binc = 1.
bpa0 = 0
cinc = 7

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
        planar[j, 5] = bmin0+(i*naxis3+j)*cinc
    gauprops.append(planar)

incubi =  ['demo_incubus1.fits', 'demo_incubus2.fits']

# Create cubes
beach.Beach(verb = False).createstcubes(gauprops = gauprops, outcubi = incubi, naxis1 = naxis1, naxis2 = naxis2, ctype3='VRAD')

# The real demo begins here, we have a list of input fitsfile names (incubi)
# and output fitsfile names (outcubi)
# The following will convolve all images to the same common beam

outcubi = ['test3_outcubus1.fits', 'test3_outcubus2.fits']
params = { 'cubes': incubi,
           'tra_fitsnames': outcubi
          }
beach.Beach(cubenames = params['cubes'], tra_fitsnames = params['tra_fitsnames'], tra_overwrite = True)
print('Created output cubes {:s} and {:s}'.format(outcubi[0], outcubi[1]))
print()
