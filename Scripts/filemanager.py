import numpy as np
import scipy as sp
from astropy.io import fits # Loading PSF with astropy module, that offers a loading function for Flexible Image Transport Image

from plotLib import saveImplot, saveImsubplots, savePlot



## Open image
from PIL import Image

def openImage(imageName: str):
    """
    Open an image file and convert it to a numpy array.

    Args:
        imageName (str): The name of the image file to open.

    Returns:
        np.ndarray: The image data as a numpy array.
        tuple: The shape of the image data.
    """
    
    image = Image.open("./data/" + imageName)
    image_data = np.array(image)
    
    ### Plot image for visual check
    saveImplot(image_data, f"Opened image: {imageName}", filename=f"opened_image_{imageName}", save=True, plot=False, logScale=False)
    shape = image_data.shape
    
    return image_data, shape

def updatePsfSize(psf, shape, verbose=False):
    """
    Update the size of the PSF to match the shape of the object.

    Args:
        psf (np.ndarray): The PSF to be updated.
        shape (tuple): The desired shape (height, width) for the PSF.
        verbose (bool, optional): If True, save and plot the updated PSF. Defaults to False.

    Returns:
        np.ndarray: The updated PSF with the desired shape.
    """
    psf = psf[150:350, 150:350] # Crop the PSF to the desired size
    res = np.pad(psf, ((2*shape[0] - psf.shape[0])//2, (2*shape[1] - psf.shape[1])//2), 'constant')
    res = sp.fft.ifftshift(res)
    res = res[:shape[0], :shape[1]]
    if verbose:
        saveImplot(res, "Updated PSF", filename="updated_psf", save=False, plot=False)
        
    return res



def createObject(nPoints = 2048, sampling_rate = 1, p=2, k=1, rho_0=0.01):
    """
    Create a synthetic object using a Gaussian random field.

    Args:
        nPoints (int, optional): Number of points in each dimension of the object. Defaults to 2048.
        sampling_rate (int, optional): Sampling rate for the frequency grid. Defaults to 1.
        p (int, optional): Power parameter for the DSP computation. Defaults to 2.
        k (int, optional): Scaling factor for the DSP computation. Defaults to 1.
        rho_0 (float, optional): Characteristic frequency for the DSP computation. Defaults to 0.01.

    Returns:
        np.ndarray: The generated synthetic object.
        np.ndarray: The DSP of the object.
    """
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

def openPsfs():
    """
    Open and read PSF data from a FITS file.

    Returns:
        np.ndarray: The PSF data read from the FITS file.
    """
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

def context2Psf(data, N=20):
    """
    Prepare the context for PSF processing.

    Args:
        data (np.ndarray): The PSF data read from the FITS file.
        N (int, optional): Number of points for the alpha axis. Defaults to 20.

    Returns:
        tuple: A tuple containing the alpha axis, Jmap axis, and the three PSFs.
    """
    psf_1, psf_2, psf_3 = data[0], data[1], data[2]
    xAxis = [round(i /(N),4) for i in range(N + 1)] # Define the alpha axis
    yAxis = np.zeros(N + 1) # Define the Jmap axis
    
    return xAxis, yAxis, psf_1, psf_2, psf_3, N 
    
def createImage(alpha_0, oObject, psf_1, psf_2, verbose=False):
    """
    Create a synthetic image by combining the object with two PSFs and adding noise.

    Args:
        alpha_0 (float): Weighting factor between psf_1 and psf_2.
        oObject (np.ndarray): The synthetic object.
        psf_1 (np.ndarray): The first PSF.
        psf_2 (np.ndarray): The second PSF.
        verbose (bool, optional): If True, print additional information. Defaults to False.

    Returns:
        np.ndarray: The synthetic image with noise.
        float: The standard deviation of the noise.
    """
    simImageNoised, sigma = createSyntheticImage(oObject, psf_1, psf_2, alpha_0, verbose=verbose)

    return simImageNoised, sigma
