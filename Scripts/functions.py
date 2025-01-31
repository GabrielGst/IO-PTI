"""
This module provides functions to compute high resolution images using estimators' theroy.

Functions:
    JmapCriterion2psf(oObject, psf_1, psf_2, alpha, verbose=False):
       Return the Jmap criterion value. The joint MAP solution is the value of alpha that minimize Jmap. 
"""


from plotLib import saveImsubplots, saveImplot
import numpy as np
import scipy as sp

def createSyntheticImage(oObject, psf_1, psf_2, alpha_0, verbose=False):
    """Create the synthetic image using the object, the two PSFs and the alpha value. The PSF is a linear combination of the two PSFs. The simulated image is computed using the object and the PSF. The noise is added to the simulated image. The noised image is returned.

    Args:
        oObject (float): Object derived from the dsp of the object.
        psf_1 (float): First PSF imported.
        psf_2 (float): Second PSF imported.
        alpha_0 (float): Between 0 and 1, the weight of the first PSF in the final PSF.
        verbose (bool, optional): Show and save the generated images (choice defined within the function). Defaults to False.

    Returns:
        tuple: (float, float) The noised image and the standard deviation of the noise.
    """
    
    # Define the synthetic image
    ## Define PSF
    psf = alpha_0 * psf_1 + (1 - alpha_0) * psf_2 # The final PSF is a linear combination of the two PSFs
    psfPadded = np.pad(psf, (768, 768), 'constant') # Object is of size 1024x1024, derived from a dsp of size 2048x2048
    
    ### Plot psf for quality check
    if verbose:
        saveImplot(psfPadded, f"PSF of the synthetic image for alpha = {alpha_0}", filename=f"PSF of the synthetic image for alpha = {alpha_0}", save=True, plot=False)


    ## Shift and crop the PSF to corresponding object size and coordinates reference
    psfPadded = sp.fft.ifftshift(psfPadded)
    psfPadded = psfPadded[:1024, :1024]
    
    
    # ## Compute the Fourier transform of the PSF
    # tfPsf = sp.fft.fft2(psfPadded) 
     
    # ### Plot psf for quality check
    # if verbose:
    #     saveImplot(np.abs(tfPsf), f"TF of the PSF for alpha = {alpha_0}", filename=f"TF of PSF for alpha = {alpha_0}", save=True, plot=True)


    ## Compute the simulated image
    simImage = sp.fft.ifft2(sp.fft.fft2(oObject) * sp.fft.fft2(psfPadded)) 
    sigma = 0.01 * np.max(oObject)
    noise = sigma * np.random.randn(*simImage.shape) # Noise corresponding roughly to photon noise for an average of 10 000 photons/pixels
    simImageNoised = simImage + noise
    
    ### Plot simulated object, noise, image for quality check
    if verbose:
        saveImsubplots([oObject, simImage.real, noise, simImageNoised.real], ["Object", "Image", "Noise", "Noised Image"], filename="simulated_object_image_noisedImage", save=True, plot=True)


    return simImageNoised, sigma

def jointMAP2psf(oObject, dspObject, psf_1, psf_2, alpha, simImageNoised, sigma):
    """Performs the joint estimation of the object and the PSF using the Jmap criterion. The criterion is computed for the given object, the two PSFs and the alpha value. The simulated image is computed using the object and the PSF. The noise is added to the simulated image. The noise and signal DSPs are computed. The Jmap criterion is computed and returned.

    Args:
        oObject (float): Object derived from the dsp of the object.
        dspObject (float): DSP of the object.
        psf_1 (float): First PSF imported.
        psf_2 (float): Second PSF imported.
        alpha (float): Between 0 and 1, the tested weight of the first PSF in the final PSF.
        simImageNoised (float): Synthetic image for which the value of alpha is searched.
        sigma (float): noise standard deviation.

    Returns:
        float: Jmap criterion value for the selected alpha.
    """
       
    ## Tested alpha and PSF
    h = alpha * psf_1 + (1 - alpha) * psf_2 # The final PSF is a linear combination of the two PSFs
    hPad = np.pad(h, (768, 768), 'constant') # Object is of size 1024x1024, derived from a dsp of size 2048x2048
    
    
    ## Shift and crop the PSF to corresponding object size and coordinates reference
    hPad = sp.fft.ifftshift(hPad)
    hPad = hPad[:1024, :1024]


    ## Compute the noise and signal DSPs
    dspObject = dspObject[:1024,:1024] # The DSP of the object is known
    N =  simImageNoised.shape[0] * simImageNoised.shape[1] # Number of pixels
    dspNoise = sigma**2 * N # dspNoise.mean() # Average noise power, the noise is white so it is constant, one could have taken dspNoise[0,0] = sigma^2 ici attention à la normalisation avec le nombre de pixel pour garder l'énergie constante
    meanObject = oObject.mean() * np.ones(oObject.shape)


    ## Compute the Jmap criterion
    Jmap = (
        0.5 * N * np.log(dspNoise) + # Here we removed the N**2 factor because it led to a very large negative offset (1e12 order of magnitude). I think it has already been taken into account in the definition of the dsp    
        0.5 * np.sum(np.log(dspObject)) +
        0.5 * np.sum(
            np.abs(sp.fft.fft2(simImageNoised) - sp.fft.fft2(hPad) * sp.fft.fft2(meanObject))**2 / # Equivalent to np.abs(sp.fft.fft2(noise))**2
            dspObject /
            (np.abs(sp.fft.fft2(hPad))**2 + dspNoise / dspObject) # Q: Is there an impact of using the padded PSF here ?
        )
    )
    
    return Jmap, dspNoise, hPad, N

def marginalML2psf(oObject, dspObject, psf_1, psf_2, alpha, simImageNoised, sigma,):
    """Performs the marginal estimation of the object and the PSF using the Jmap criterion. The criterion is computed for the given object, the two PSFs and the alpha value. The simulated image is computed using the object and the PSF. The noise is added to the simulated image. The noise and signal DSPs are computed. The Jmap criterion is computed and returned.

    Args:
        oObject (float): Object derived from the dsp of the object.
        dspObject (float): DSP of the object.
        psf_1 (float): First PSF imported.
        psf_2 (float): Second PSF imported.
        alpha (float): Between 0 and 1, the tested weight of the first PSF in the final PSF.
        simImageNoised (float): Synthetic image for which the value of alpha is searched.
        sigma (float): noise standard deviation.

    Returns:
        (float, float, int): (Jmap, dspNoise, N) Jmap criterion value for the selected alpha, noise DSP and number of pixels.
    """
    
    ## Retrieve Jmap criterion and other parameters
    Jmap, dspNoise, hPad, N = jointMAP2psf(oObject, dspObject, psf_1, psf_2, alpha, simImageNoised, sigma)
    
    
    ## Crop the dspObject to the same size as the object
    dspObject = dspObject[:1024,:1024]


    ## Compute the Jml criterion
    Jml = (
        Jmap +
        0.5 * np.sum(np.log(np.abs(sp.fft.fft2(hPad))**2 + dspNoise / dspObject)) - 
        0.5 * N * np.log(dspNoise) # Here we removed the N**2 factor because it led to a very large negative offset (1e12 order of magnitude). I think it has already been taken into account in the definition of the dsp
    )
    
    return Jml, dspNoise, N

def estimatedImage(simImageNoised, oObject, dspObject, dspNoise, alpha, psf_1, psf_2):
    """Reconstruct an estimated image based on the alpha resulting from the Marginal estimation. The estimated image is computed using the object, the simulated image, the noise and the two PSFs. The final PSF is a linear combination of the two PSFs.

    Args:
        simImageNoised (float): Synthetic image for which the value of alpha is searched.
        oObject (float): Object derived from the dsp of the object.
        dspObject (float): DSP of the object.
        dspNoise (float): DSP of the noise.
        alpha (float): Estimated value of alpha.
        psf_1 (float): First PSF imported.
        psf_2 (float): Second PSF imported.
    """
    
    ## Define the estimated PSF and average of the object
    psf = alpha * psf_1 + (1 - alpha) * psf_2 # The final PSF is a linear combination of the two PSFs
    oM = np.mean(oObject) * np.ones(oObject.shape)
    
    ## Compute the estimated image
    im = (
        (sp.fft.fft2(psf.conj()) * sp.fft.fft2(simImageNoised) + dspNoise / dspObject * oM) /
         (np.abs(sp.fft.fft2(psf))**2 + dspNoise / dspObject)
    )
    
    ### Plot psf for quality check
    saveImplot(im, f"Estimated object for estimated alpha = {alpha}", filename=f"Estimated object for estimated alpha = {alpha}", save=True, plot=True)

