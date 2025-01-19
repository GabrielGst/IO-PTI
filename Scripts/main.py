import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from plotLib import saveImplot, saveImsubplots, savePlot

# Define the path to save figures in plotLib SAVE_PATH = 'IO-PTI/outputs/' # Q: Why do we need to use the absolute path ?





# Define the gaussian field in polar coordinates (freq.)
## Define the polar grid
nPoints = 2048 # Number of points
sampling_rate = 1  # Sampling rate
freq_vector = sp.fft.fftfreq(nPoints, d=1/sampling_rate) # Frequency vector for FFT, Numpy array
xFreq_array, yFreq_array = np.meshgrid(freq_vector, freq_vector) # Generating 2D array of frequencies from freq_vector


## Convert Cartesian (X, Y) to Polar (r, theta)
rho = np.hypot(xFreq_array, yFreq_array)  # Radial distance
theta = np.arctan2(xFreq_array, yFreq_array)  # Angle in radians


### Plot modulus and phase for quality check
saveImplot(rho, "Polar modulus", filename="polar_modulus", save=True, plot=True)
saveImplot(theta, "Polar phase", filename="polar_phase", save=True, plot=True)


## Computes the dsp of the object
p = 2
k = 1
rho_0 = 0.01

dspObject = k /(1 + pow(rho/rho_0,p))
profile = dspObject[0,:512] # Profile of the radial distance
freqProfile = freq_vector[:512] # Profile of the frequency vector


### Plot dspObject for quality check
savePlot(freqProfile, profile, "Profile of the radial distance", filename="profile_radial_distance", save=True, plot=True, logX=True, logY=True)
saveImplot(dspObject, f"dspObject for rho_0 = {rho_0}", filename=f"dspObject for rho_0 = {rho_0}", save=True, plot=True, logScale=True)


## Generate the random variable
tfObject = np.sqrt(dspObject) * np.random.randn(nPoints, nPoints) # Equivalently : np.sqrt(dspObject) * np.random.randn(nPoints, nPoints)



## Derive the object
oObject = np.real(sp.fft.ifft2(tfObject))
oShapes = oObject.shape
oObject = oObject[:int(oShapes[0]/2),:int(oShapes[1]/2)] # Select first quadrant, int is equivalent to floor
oObject = oObject - np.min(oObject) # Remove the offset


### Plot oObject for quality check
saveImplot(oObject, f"Object for rho_0 = {rho_0}", filename=f"Object for rho_0 = {rho_0}", save=True, plot=True, logScale=False)



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
saveImsubplots(data, ["PSF 1", "PSF 2", "PSF 3"], filename="psfs", save=True, plot=True)





# Bayesian estimation : Derive alpha from 2 extreme PSF and criterion
from functions import JmapCriterion2psf
from tqdm.auto import tqdm

psf_1, psf_2, psf_3 = data[0], data[1], data[2]
N = 20
xAxis = [i /(N + 1) for i in range(N + 1)] # Define the alpha axis
yAxis = np.zeros(N + 1) # Define the Jmap axis

pbar = tqdm(desc="Jmap computation", total=100, unit_scale=True)

for i, alpha in enumerate(xAxis):
    yAxis[i] = JmapCriterion2psf(oObject, dspObject, psf_1, psf_3, alpha, verbose=True if i==1 else False)
    pbar.update(round(100 /(N + 1)))
pbar.close()

### Plot Jmap for quality check
savePlot(xAxis, yAxis, "Jmap criterion", filename="Jmap_criterion", save=True, plot=True, logX=False, logY=False)