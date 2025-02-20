from tqdm.auto import tqdm
import numpy as np
import scipy as sp
from plotLib import saveImplot, saveImsubplots, savePlot

## Joint estimation
from functions import jointMAP2psf

def jointEstimation(oObject, dspObject, psf_1, psf_3, simImageNoised, sigma, xAxis, yAxis, N, name):
    """
    Perform joint estimation of the PSF using the MAP criterion.

    Args:
        oObject (np.ndarray): The synthetic object.
        dspObject (np.ndarray): The DSP of the object.
        psf_1 (np.ndarray): The first PSF.
        psf_3 (np.ndarray): The third PSF.
        simImageNoised (np.ndarray): The synthetic image with noise.
        sigma (float): The standard deviation of the noise.
        xAxis (list): The alpha axis.
        yAxis (np.ndarray): The Jmap axis.
        N (int): Number of points for the alpha axis.
    """
    pbar = tqdm(desc="Jmap computation", total=N, unit_scale=True)

    for i, alpha in enumerate(xAxis):
        res = jointMAP2psf(oObject, dspObject, psf_1, psf_3, alpha, simImageNoised, sigma)
        yAxis[i] = res[0]
        # print(f"Jmap Value for i = {i}: {yAxis[i]}")
        pbar.update() #round(100 /(N - 1))
    pbar.close()

    ### Plot Jmap for quality check
    savePlot(xAxis, yAxis, "MAP-estimated alpha", filename=f"MAP-estimated alpha {name}", save=True, plot=True, logX=False, logY=False)


## Marginal estimation
from functions import marginalML2psf

def marginalEstimation(oObject, dspObject, psf_1, psf_3, simImageNoised, sigma, xAxis, yAxis, N, name):
    """
    Perform marginal estimation of the PSF using the ML criterion.

    Args:
        oObject (np.ndarray): The synthetic object.
        dspObject (np.ndarray): The DSP of the object.
        psf_1 (np.ndarray): The first PSF.
        psf_3 (np.ndarray): The third PSF.
        simImageNoised (np.ndarray): The synthetic image with noise.
        sigma (float): The standard deviation of the noise.
        xAxis (list): The alpha axis.
        yAxis (np.ndarray): The Jml axis.
        N (int): Number of points for the alpha axis.
    """
       
    pbar = tqdm(desc="ML-estimated alpha", total=N, unit_scale=True)

    for i, alpha in enumerate(xAxis):
        res = marginalML2psf(oObject, dspObject, psf_1, psf_3, alpha, simImageNoised, sigma)
        yAxis[i] = res[0]
        # print(f"Jml Value for i = {i}: {yAxis[i]}")
        pbar.update() # round(100 /(N - 1))
    pbar.close()

    ### Plot Jmap for quality check
    savePlot(xAxis, yAxis, f"ML-estimated alpha {name}", filename=f"ML-estimated alpha {name}", save=True, plot=True, logX=False, logY=False)

    alpha_min = xAxis[yAxis.argmin()]
    print(f"Minimal argument of estimation function : {alpha_min}")
    
    return alpha_min

## Estimation

def estimatedObject(image, oObject, dspObject, dspNoise, alpha, psf_1, psf_2, name):
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
    o = sp.fft.ifft2(
        (sp.fft.fft2(psf).conj() * sp.fft.fft2(image) + dspNoise / dspObject * sp.fft.fft2(oM)) /
        (np.abs(sp.fft.fft2(psf))**2 + dspNoise / dspObject)
    )
    
    ### Plot psf for quality check
    o = o.real
    o = o - o.min()
    # o[o > 1e-3] = 0
    # o = o.real - o.real.mean()
    saveImplot(o, f"ML-estimated object ({name}) for estimated alpha = {alpha}", filename=f"ML-estimated object ({name}) for estimated alpha = {alpha}", save=True, plot=True, logScale=False)
    saveImsubplots([oObject.real, o], [f"Original ({name})", f"Estimated object ({name})"], filename=f"original vs estimated object ({name})", save=True, plot=True)
    return o
