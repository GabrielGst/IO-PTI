import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from astropy.io import fits # Loading PSF with astropy module, that offers a loading function for Flexible Image Transport Image
from tqdm.auto import tqdm


from plotLib import saveImplot, saveImsubplots, savePlot

# Define the path to save figures in plotLib SAVE_PATH = 'IO-PTI/outputs/' # Q: Why do we need to use the absolute path ?

def createObject(nPoints = 2048, sampling_rate = 1, p=2, k=1, rho_0=0.01):
    # Define the gaussian field in polar coordinates (freq.)
    ## Define the polar grid
    freq_vector = sp.fft.fftfreq(nPoints, d=1/sampling_rate) # Frequency vector for FFT, Numpy array
    xFreq_array, yFreq_array = np.meshgrid(freq_vector, freq_vector) # Generating 2D array of frequencies from freq_vector


    ## Convert Cartesian (X, Y) to Polar (r, theta)
    rho = np.hypot(xFreq_array, yFreq_array)  # Radial distance
    theta = np.arctan2(xFreq_array, yFreq_array)  # Angle in radians


    ### Plot modulus and phase for quality check
    saveImplot(rho, "Polar modulus", filename="polar_modulus", save=True, plot=False)
    saveImplot(theta, "Polar phase", filename="polar_phase", save=True, plot=False)


    ## Computes the dsp of the object
    dspObject = k /(1 + pow(rho/rho_0,p))
    profile = dspObject[0,:1024] # Profile of the radial distance
    freqProfile = freq_vector[:1024] # Profile of the frequency vector


    ### Plot dspObject for quality check
    savePlot(freqProfile, profile, "Profile of the radial distance", filename="profile_radial_distance", save=True, plot=False, logX=True, logY=True)
    saveImplot(dspObject, f"dspObject for rho_0 = {rho_0}", filename=f"dspObject for rho_0 = {rho_0}", save=True, plot=False, logScale=True)


    ## Generate the random variable
    tfObject = np.sqrt(dspObject) * np.random.randn(nPoints, nPoints) # Equivalently : np.sqrt(dspObject) * np.random.randn(nPoints, nPoints)



    ## Derive the object
    oObject = np.real(sp.fft.ifft2(tfObject))
    oShapes = oObject.shape
    oObject = oObject[:int(oShapes[0]/2),:int(oShapes[1]/2)] # Select first quadrant, int is equivalent to floor
    oObject = oObject - np.min(oObject) # Remove the offset


    ### Plot oObject for quality check
    saveImplot(oObject, f"Object for rho_0 = {rho_0}", filename=f"Object for rho_0 = {rho_0}", save=True, plot=False, logScale=False)

    return oObject, dspObject

def openPsfs(system: str):
    if system == "windows":
        filename = 'IO-PTI/data/3psfs_zeroPiSur2EtPi.fits' # Open the FITS file on windows
    elif system == "ubuntu":
        filename = './data/3psfs_zeroPiSur2EtPi.fits' # Open the FITS file on ubuntu

    with fits.open(filename) as hdul:
        # hdul is a list-like object containing all HDUs (Header Data Units)
        hdul.info()  # Display information about the HDUs
        primary_hdu = hdul[0] # Access the primary HDU (header and data)
        header = primary_hdu.header # The header
        # print("Header:", header)
        data = primary_hdu.data # The data (typically image data or a 2D array)
        # print("Data:", data)
    
    ### Plot psfs for visual check
    saveImsubplots(data, ["PSF 1", "PSF 2", "PSF 3"], filename="psfs", save=True, plot=False)


    return data


# Bayesian estimation : Derive alpha from 2 extreme PSF and criterion
from functions import createSyntheticImage

def createImage(alpha, data, oObject):
    psf_1, psf_2, psf_3 = data[0], data[1], data[2]
    N = 20
    xAxis = [round(i /(N + 1),4) for i in range(N + 1)] # Define the alpha axis
    yAxis = np.zeros(N + 1) # Define the Jmap axis

    alpha_0 = alpha
    simImageNoised, sigma = createSyntheticImage(oObject, psf_1, psf_2, alpha_0, verbose=True)

    return simImageNoised, sigma, xAxis, yAxis, psf_1, psf_3, N


## Joint estimation
from functions import jointMAP2psf

def jointEstimation(oObject, dspObject, psf_1, psf_3, simImageNoised, sigma, xAxis, yAxis, N):    
    pbar = tqdm(desc="Jmap computation", total=100, unit_scale=True)

    for i, alpha in enumerate(xAxis):
        res = jointMAP2psf(oObject, dspObject, psf_1, psf_3, alpha, simImageNoised, sigma)
        yAxis[i] = res[0]
        # print(f"Jmap Value for i = {i}: {yAxis[i]}")
        pbar.update(round(100 /(N + 1)))
    pbar.close()

    ### Plot Jmap for quality check
    savePlot(xAxis, yAxis, "Jmap criterion", filename="Jmap_criterion", save=True, plot=True, logX=False, logY=False)


## Marginal estimation
from functions import marginalML2psf

def marginalEstimation(oObject, dspObject, psf_1, psf_3, simImageNoised, sigma, xAxis, yAxis, N):    
    pbar = tqdm(desc="Jml computation", total=100, unit_scale=True)

    for i, alpha in enumerate(xAxis):
        res = marginalML2psf(oObject, dspObject, psf_1, psf_3, alpha, simImageNoised, sigma)
        yAxis[i] = res[0]
        # print(f"Jml Value for i = {i}: {yAxis[i]}")
        pbar.update(round(100 /(N + 1)))
    pbar.close()

    ### Plot Jmap for quality check
    savePlot(xAxis, yAxis, "Jml criterion", filename="Jml_criterion", save=True, plot=True, logX=False, logY=False)


def main():
    oObject, dspObject = createObject()
    psfs = openPsfs("ubuntu")
    simImageNoised, sigma, xAxis, yAxis, psf_1, psf_3, N = createImage(0.3, psfs, oObject)
    jointEstimation(oObject, dspObject, psf_1, psf_3, simImageNoised, sigma, xAxis, yAxis, N)
    marginalEstimation(oObject, dspObject, psf_1, psf_3, simImageNoised, sigma, xAxis, yAxis, N)

    
if __name__ == "__main__":
    main()
