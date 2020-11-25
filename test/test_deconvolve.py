from equolver.beach import Beach
import astropy.io as astropy_io
from astropy.io import fits
import numpy as np
from numpy.testing import assert_almost_equal
import pyfftw
import pytest

# Under construction


@pytest.fixture
def cube():
    return np.ones((100, 100, 128))

@pytest.fixture(params=[16, 32])
def freq(request):
    return np.linspace(.856e9, 2*.856e9, request.param)


@pytest.mark.parametrize("corrs", [
    ("XX", "XY", "YX", "YY"),
    ("RR", "RL", "LR", "LL"),
])
@pytest.mark.parametrize("freq", [8, 64], indirect=True)
#@pytest.mark.parametrize("freq", [np.linspace(.856e9, 2*.856e9, 16)])
def test_example(freq, corrs, cube):
    print(corrs, freq.size, cube.sum())


@pytest.mark.parametrize("point_source", ["point_source.fits"])
@pytest.mark.parametrize("gaussian_at_centre", ["gaussian_at_centre.fits"])
@pytest.mark.parametrize("gaussian_at_origin", ["gaussian_at_origin.fits"])
@pytest.mark.parametrize("real_fft_conv", ["real_fft_conv.fits"])
@pytest.mark.parametrize("real_fft_conv_calc", ["real_fft_conv_calc.fits"])
@pytest.mark.parametrize("reconvolve_input_image", ["reconvolve_input_image.fits"])
@pytest.mark.parametrize("reconvolve_output_image", ["reconvolve_output_image.fits"])
def test_convoltests(point_source, gaussian_at_centre, gaussian_at_origin, real_fft_conv, real_fft_conv_calc, reconvolve_input_image, reconvolve_output_image):
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

    beach = Beach()


    hdu[0].data[:] = beach._gaussian_2dp( naxis1 =
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

    hdu[0].data[:] = beach._gaussian_2dp( naxis1 =
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

    kernel = beach._gaussian_2dp(naxis1 = target.shape[1], naxis2 =
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

    ikernel = beach._igaussian_2dp(naxis1 = target.shape[1], naxis2
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

    hdu[0].data[:] = beach._gaussian_2dp( naxis1 =
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


    orig_centre = hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2]
    orig_minor_half_power = hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)]
    orig_major_half_power = hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)]

    print('Original result at centre', orig_centre)
    print('Original result at half power', orig_minor_half_power)
    print('Original result at half power', orig_major_half_power)

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

    #ikernel = beach._igaussian_2dp(naxis1 = target.shape[1], naxis2
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
    ikernel = beach._igaussian_2dp(naxis1 = target.shape[1], naxis2
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

    final_centre = hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2]
    final_minor_half_power = hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)]
    final_major_half_power = hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)]

    assert_almost_equal(orig_centre, final_centre)
    assert_almost_equal(orig_minor_half_power, final_minor_half_power)
    assert_almost_equal(orig_major_half_power, final_major_half_power)

    print('Final result at centre', final_centre)
    print('Final result at half power', final_minor_half_power)
    print('Final result at half power', final_major_half_power)
    hdu[0].data[hdu[0].data.shape[0]//2, hdu[0].data.shape[1]//2] += 2
    hdu[0].data[hdu[0].data.shape[0]//2+int(ymin),hdu[0].data.shape[1]//2+int(xmin)] += 2.
    hdu[0].data[hdu[0].data.shape[0]//2+int(ymaj),hdu[0].data.shape[1]//2+int(xmaj)] += 2.
    hdu.writeto(reconvolve_output_image, overwrite = True)
    print('Image to be found at', reconvolve_output_image)
    hdu.close()
    print('')

if __name__ == '__main__':
    convoltests()
