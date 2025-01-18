import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# Define the gaussian field in polar coordinates (freq.)
## Define the polar grid
nPoints = 2048 # Number of points
sampling_rate = 1  # Sampling rate
freq_vector = sp.fft.fftfreq(nPoints, d=1/sampling_rate) # Frequency vector for FFT, Numpy array
freq_vector = sp.fft.fftshift(freq_vector) # perform a fftshift
xFreq_array, yFreq_array = np.meshgrid(freq_vector, freq_vector) # Generating 2D array of frequencies from freq_vector

## Computes the distance and plot it for quality check. Frequencies increases from the center of the figure
# distCartesian = np.sqrt(xFreq_array**2 + yFreq_array**2)

# titlePlot = "Cartesian distance"
# plt.figure()
# plt.imshow(distCartesian, cmap='gray')
# plt.title(titlePlot)
# plt.colorbar()
# plt.savefig(f'{titlePlot}.png') # Save the figure for remote shell
# # plt.show()

# Convert Cartesian (X, Y) to Polar (r, theta)
rho = np.hypot(xFreq_array, yFreq_array)  # Radial distance
theta = np.arctan2(xFreq_array, yFreq_array)  # Angle in radians


## Plot modulus and phase for quality check
# titlePlot = "Polar modulus"
# plt.figure()
# plt.imshow(rho, cmap='gray')
# plt.title(titlePlot)
# plt.colorbar()
# plt.savefig(f'{titlePlot}.png') # Save the figure for remote shell
# # plt.show()

# titlePlot = "Polar phase"
# plt.figure()
# plt.imshow(theta, cmap='gray')
# plt.title(titlePlot)
# plt.colorbar()
# plt.savefig(f'{titlePlot}.png') # Save the figure for remote shell
# # plt.show()


## Computes the dsp of the object
p = 2
k = 1
rho_0 = 0.01

dspObject = k /(1+pow(rho/rho_0,p))

## Plot dspObject for quality check
from matplotlib.colors import LogNorm

# titlePlot = f"dspObject for rho_0 = {rho_0}"
# plt.figure()
# plt.imshow(dspObject, cmap='gray', norm=LogNorm(vmin=0.01, vmax=1))
# plt.title(titlePlot)
# plt.colorbar()
# plt.savefig(f'{titlePlot}.png') # Save the figure for remote shell
# # plt.show()


## Generate the random variable
tfObject = np.sqrt(dspObject) * np.random.randn(nPoints, nPoints) # Equivalently : np.sqrt(dspObject) * np.random.randn(nPoints, nPoints)

## Derive the object
oObject = np.real(sp.fft.ifft2(sp.fft.ifftshift(tfObject)))
oShapes = oObject.shape
oObject = oObject[:int(oShapes[0]/2),:int(oShapes[1]/2)] # Select first quadrant, int is equivalent to floor
oObject = oObject - np.min(oObject) # Remove the offset

## Plot oObject for quality check
# titlePlot = f"Object for rho_0 = {rho_0}"
# plt.figure()
# plt.imshow(oObject, cmap='gray')
# plt.title(titlePlot)
# plt.colorbar()
# plt.savefig(f'{titlePlot}.png') # Save the figure for remote shell
# # plt.show()


# Loading PSF with astropy module, that offers a loading function for Flexible Image Transport Image
from astropy.io import fits

filename = '3psfs_zeroPiSur2EtPi.fits' # Open the FITS file

with fits.open(filename) as hdul:
     # hdul is a list-like object containing all HDUs (Header Data Units)
    hdul.info()  # Display information about the HDUs

    # Access the primary HDU (header and data)
    primary_hdu = hdul[0]

    # The header
    header = primary_hdu.header
    # print("Header:", header)

    # The data (typically image data or a 2D array)
    data = primary_hdu.data
    # print("Data:", data)

# titlePlot = f"PSFs"

# plt.figure(figsize=(12, 4)) # Create a 1x3 grid for the subplots

# # First subplot (1 row, 3 columns, first plot)
# plt.subplot(1, 3, 1)  # (rows, columns, index)
# plt.imshow(data[0], cmap='gray')
# plt.title("file 1")

# # Second subplot (1 row, 3 columns, second plot)
# plt.subplot(1, 3, 2)
# plt.imshow(data[1], cmap='gray')
# plt.title("file 2")

# # Third subplot (1 row, 3 columns, third plot)
# plt.subplot(1, 3, 3)
# plt.imshow(data[2], cmap='gray')
# plt.title("file 3")

# # Adjust layout for better spacing
# plt.tight_layout()

# # Show the plots
# plt.savefig(f'{titlePlot}.png') # Save the figure for remote shell
# # plt.show()


# Bayesian estimation : Derive alpha from 2 extreme PSF and criterion
psf_1, psf_2, psf_3 = data[0], data[1], data[2]

## Define PSF
alpha = 0.3
psf = alpha * psf_1 + (1 - alpha) * psf_3

## Plot psf for quality check
titlePlot = f"Resulting psf for alpha = {alpha}"
plt.figure()
plt.imshow(psf, cmap='gray')
plt.title(titlePlot)
plt.colorbar()
plt.savefig(f'{titlePlot}.png') # Save the figure for remote shell
# plt.show()