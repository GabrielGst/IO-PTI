"""
This module provides functions to compute high resolution images using estimators' theroy.

Functions:
    JmapCriterion2psf(oObject, psf_1, psf_2, alpha, verbose=False):
       Return the Jmap criterion value. The joint MAP solution is the value of alpha that minimize Jmap. 
"""


from plotLib import saveImsubplots, saveImplot
import numpy as np
import scipy as sp

def jointMAP2psf(oObject, dspObject, psf_1, psf_2, alpha, verbose=False):
    """Performs the joint estimation of the object and the PSF using the Jmap criterion. The criterion is computed for the given object, the two PSFs and the alpha value. The simulated image is computed using the object and the PSF. The noise is added to the simulated image. The noise and signal DSPs are computed. The Jmap criterion is computed and returned.

    Args:
        oObject (float): Array representing the object to be estimated.
        psf_1 (float): Array representing the first PSF.
        psf_2 (float): Array representing the second PSF.
        alpha (float): Scalar representing the weight of the first PSF in the final PSF.
        verbose (bool, optional): Display and save intermediate plots. Defaults to False.

    Returns:
        float: Jmap criterion value. The joint MAP solution is the value of alpha that minimize Jmap.
    """
    ## Define PSF
    # alpha = 0.3
    psf = alpha * psf_1 + (1 - alpha) * psf_2


    ### Plot psf for quality check
    if verbose:
        saveImplot(psf, f"Resulting psf for alpha = {alpha}", filename=f"Resulting psf for alpha = {alpha}", save=False, plot=True)



    ## Joint estimation
    # from scipy.signal import convolve2d

    psfPadded = np.pad(psf, (256, 256), 'wrap') # Use the padded psf for the Fourier transform
    psfPadded = sp.fft.ifftshift(psfPadded)
    simImage = sp.fft.ifft2(sp.fft.fft2(oObject) * sp.fft.fft2(psfPadded))# Q: Should the psf be fftshifted here ?
    # simImage[:512, :512]
    sigma = 0.01 * np.max(oObject)
    noise = sigma * np.random.randn(*simImage.shape) # Noise corresponding roughly to photon noise for an average of 10 000 photons/pixels
    simImageNoised = simImage + noise


    ### Plot simulated object, noise, image for quality check
    if verbose:
        saveImsubplots([oObject, simImage.real, noise, simImageNoised.real], ["Object", "Image", "Noise", "Noised Image"], filename="simulated_object_image_noisedImage", save=True, plot=True)



    ## Compute the noise and signal DSPs
    # dspNoise = np.abs(sp.fft.ifft2(sp.fft.fft2(noise) * sp.fft.fft2(noise).conj()))**2 # Correlation therorem
    # dspObject = sp.fft.ifft2(sp.fft.fft2(oObject) * sp.fft.fft2(oObject).conj()).real # Correlation therorem
    dspObject = dspObject[:1024,:1024] # The DSP of the object is known
    

    ### Plot simulated object, noise, image for quality check
    # if verbose:
    #     saveImsubplots([dspNoise, dspObject], ["DSP of the synthetic noise", "DSP of the synthetic object"], filename="dsp_noise_and_object", columnNumber=2, save=True, plot=True, logScale=True)


    ## Compute the Jjmap criterion and plot it
    N =  simImage.shape[0] * simImage.shape[1] # Number of pixels
    dspNoise = sigma**2 * simImage.shape[0]**2 # dspNoise.mean() # Average noise power, the noise is white so it is constant, one could have taken dspNoise[0,0] = sigma^2 ici attention à la normalisation avec le nombre de pixel pour garder l'énergie constante

    Jmap = (
        0.5 * N * np.log(dspNoise) + # Here we removed the N**2 factor because it led to a very large negative offset (1e12 order of magnitude). I think it has already been taken into account in the definition of the dsp    
        0.5 * np.sum(np.log(dspObject)) +
        0.5 * np.sum(
            np.abs(sp.fft.fft2(simImageNoised) - sp.fft.fft2(simImage))**2 / # Equivalent to np.abs(sp.fft.fft2(noise))**2
            dspObject /
            (np.abs(sp.fft.fft2(psfPadded))**2 + dspNoise * np.reciprocal(dspObject)) # Q: Is there an impact of using the padded PSF here ?
        )
    )
    
    return Jmap, dspNoise, psfPadded, N


def marginalML2psf(oObject, dspObject, psf_1, psf_2, alpha, verbose=False):
    """Performs the joint estimation of the object and the PSF using the Jmap criterion. The criterion is computed for the given object, the two PSFs and the alpha value. The simulated image is computed using the object and the PSF. The noise is added to the simulated image. The noise and signal DSPs are computed. The Jmap criterion is computed and returned.

    Args:
        oObject (float): Array representing the object to be estimated.
        psf_1 (float): Array representing the first PSF.
        psf_2 (float): Array representing the second PSF.
        alpha (float): Scalar representing the weight of the first PSF in the final PSF.
        verbose (bool, optional): Display and save intermediate plots. Defaults to False.

    Returns:
        float: Jmap criterion value. The joint MAP solution is the value of alpha that minimize Jmap.
    """
    Jmap, dspNoise, psfPadded, N = jointMAP2psf(oObject, dspObject, psf_1, psf_2, alpha, verbose)
    dspObject = dspObject[:1024,:1024]

    Jml = (
        # Jmap +
        0.5 * np.sum(np.log(np.abs(sp.fft.fft2(psfPadded))**2 + dspNoise * np.reciprocal(dspObject))) - 
        0.5 * N * np.log(dspNoise) # Here we removed the N**2 factor because it led to a very large negative offset (1e12 order of magnitude). I think it has already been taken into account in the definition of the dsp
    )
    
    return Jml