"""
This module provides functions to compute high resolution images using estimators' theroy.

Functions:
    JmapCriterion2psf(oObject, psf_1, psf_2, alpha, verbose=False):
       Return the Jmap criterion value. The joint MAP solution is the value of alpha that minimize Jmap. 
"""


from plotLib import saveImsubplots, saveImplot
import numpy as np
import scipy as sp

def JmapCriterion2psf(oObject, psf_1, psf_2, alpha, verbose=False):
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
    simImage = sp.fft.ifft2(sp.fft.fft2(oObject) * sp.fft.fft2(psfPadded)).real # Q: Should the psf be fftshifted here ?
    simImage[:256, :256]
    noise = 0.01 * np.max(oObject) * np.random.randn(*simImage.shape) # Noise corresponding roughly to photon noise for an average of 10 000 photons/pixels
    simImageNoised = simImage + noise


    ### Plot simulated object, noise, image for quality check
    if verbose:
        saveImsubplots([oObject, simImage, noise, simImageNoised], ["Object", "Image", "Noise", "Noised Image"], filename="simulated_object_image_noisedImage", save=True, plot=True)



    ## Compute the noise and signal DSPs
    dspNoise = sp.fft.ifft2(sp.fft.fft2(noise) * sp.fft.fft2(noise).conj()).real # Correlation therorem
    dspObject = sp.fft.ifft2(sp.fft.fft2(oObject) * sp.fft.fft2(oObject).conj()).real # Correlation therorem



    ## Compute the Jjmap criterion and plot it
    N =  simImage.shape[0] * simImage.shape[1] # Number of pixels
    dspNoise = dspNoise.mean() # Average noise power, the noise is white so it is constant, one could have taken dspNoise[0,0]

    Jmap = (
        0.5 * N**2 * np.log(dspNoise) +
        0.5 * np.sum(np.log(dspObject)) +
        0.5 * np.sum(
            np.abs(sp.fft.fft2(simImageNoised) - sp.fft.fft2(simImage))**2 /
            dspObject /
            (np.abs(sp.fft.fft2(psfPadded))**2 + dspNoise / dspObject) # Q: Is there an impact of using the padded PSF here ?
        )
    )
    
    return Jmap