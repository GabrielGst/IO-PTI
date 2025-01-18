import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from plotLib import saveImplot, saveImsubplots, savePlot

# Define the path to save figures
SAVE_PATH = 'IO-PTI/outputs/' # Q: Why do we need to use the absolute path ?





# Define the gaussian field in polar coordinates (freq.)
## Define the polar grid
nPoints = 2048 # Number of points
sampling_rate = 1  # Sampling rate
freq_vector = sp.fft.fftfreq(nPoints, d=1/sampling_rate) # Frequency vector for FFT, Numpy array
freq_vector = sp.fft.fftshift(freq_vector) # perform a fftshift
xFreq_array, yFreq_array = np.meshgrid(freq_vector, freq_vector) # Generating 2D array of frequencies from freq_vector


### Computes the distance and plot it for quality check. Frequencies increases from the center of the figure
# distCartesian = np.sqrt(xFreq_array**2 + yFreq_array**2)
# saveImplot(distCartesian, "Cartesian distance", filename="cartesian_distance", save=False, plot=True)



## Convert Cartesian (X, Y) to Polar (r, theta)
rho = np.hypot(xFreq_array, yFreq_array)  # Radial distance
theta = np.arctan2(xFreq_array, yFreq_array)  # Angle in radians


### Plot modulus and phase for quality check
# saveImplot(rho, "Polar modulus", filename="polar_modulus", save=False, plot=True)
# saveImplot(theta, "Polar phase", filename="polar_phase", save=False, plot=True)


## Computes the dsp of the object
p = 2
k = 1
rho_0 = 0.01

dspObject = k /(1 + pow(rho/rho_0,p))
profile = dspObject[1024,512:] # Profile of the radial distance
freqProfile = freq_vector[512:] # Profile of the frequency vector


### Plot dspObject for quality check
# savePlot(freqProfile, profile, "Profile of the radial distance", filename="profile_radial_distance", save=False, plot=True, logX=True, logY=True)
# saveImplot(dspObject, f"dspObject for rho_0 = {rho_0}", filename=f"dspObject for rho_0 = {rho_0}", save=False, plot=True, logScale=True)


## Generate the random variable
tfObject = np.sqrt(dspObject) * np.random.randn(nPoints, nPoints) # Equivalently : np.sqrt(dspObject) * np.random.randn(nPoints, nPoints)



## Derive the object
oObject = np.real(sp.fft.ifft2(sp.fft.ifftshift(tfObject)))
oShapes = oObject.shape
oObject = oObject[:int(oShapes[0]/2),:int(oShapes[1]/2)] # Select first quadrant, int is equivalent to floor
oObject = oObject - np.min(oObject) # Remove the offset


### Plot oObject for quality check
# saveImplot(oObject, f"Object for rho_0 = {rho_0}", filename=f"Object for rho_0 = {rho_0}", save=False, plot=True, logScale=False)



## Loading PSF with astropy module, that offers a loading function for Flexible Image Transport Image
from astropy.io import fits

filename = 'IO-PTI/data/3psfs_zeroPiSur2EtPi.fits' # Open the FITS file

with fits.open(filename) as hdul:
    # hdul is a list-like object containing all HDUs (Header Data Units)
    hdul.info()  # Display information about the HDUs
    primary_hdu = hdul[0] # Access the primary HDU (header and data)
    header = primary_hdu.header # The header
    # print("Header:", header)
    data = primary_hdu.data # The data (typically image data or a 2D array)
    # print("Data:", data)

### Plot psfs for visual check
# saveImsubplots(data, ["PSF 1", "PSF 2", "PSF 3"], filename="psfs", save=False, plot=True)





# Bayesian estimation : Derive alpha from 2 extreme PSF and criterion
psf_1, psf_2, psf_3 = data[0], data[1], data[2]



## Define PSF
alpha = 0.3
psf = alpha * psf_1 + (1 - alpha) * psf_3



### Plot psf for quality check
saveImplot(psf, f"Resulting psf for alpha = {alpha}", filename=f"Resulting psf for alpha = {alpha}", save=False, plot=True)



## Joint estimation
# from scipy.signal import convolve2d

simulatedObject = oObject
psfPadded = np.pad(psf, (256, 256), 'wrap') # Use the padded psf for the Fourier transform
simulatedImage = sp.fft.ifft2(sp.fft.fft2(simulatedObject) * sp.fft.fft2(psfPadded)).real # Q: Should the psf be fftshifted here ?
simulatedImage[:256, :256]
noise = 0.01 * np.max(simulatedObject) * np.random.randn(*simulatedImage.shape) # Noise corresponding rooughly to photon noise for an average of 10 000 photons/pixels
simulatedImageNoised = simulatedImage + noise



### Plot simulated object, noise, image for quality check
saveImsubplots([simulatedObject, simulatedImage, noise, simulatedImageNoised], ["Object", "Image", "Noise", "Noised Image"], filename="simulated_object_image_noisedImage", save=True, plot=True)


## Compute the noise and signal DSPs
noiseDsp



## Compute the Jjmap criterion and plot it


