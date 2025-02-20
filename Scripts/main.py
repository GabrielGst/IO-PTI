import numpy as np
import scipy as sp

from filemanager import updatePsfSize, openImage, openPsfs, createImage, createObject, context2Psf
from estimate import jointEstimation, marginalEstimation, estimatedObject
from plotLib import saveImplot, saveImsubplots, savePlot

# Define the path to save figures in plotLib SAVE_PATH = 'IO-PTI/outputs/' # Q: Why do we need to use the absolute path ?



def main():
    psfs = openPsfs()
    xAxis, yAxis, psf_1, psf_2, psf_3, N = context2Psf(psfs, N=40)
    
    # Synthetic object / image estimation, with two PSFs
    oObject, dspObject = createObject("Synthetic", nPoints=2048, sampling_rate=1, p=2.2, k=1, rho_0=0.01)
    psf_A, psf_C = updatePsfSize(psf_1, oObject.shape, name="Synthetic A", verbose=True), updatePsfSize(psf_3, oObject.shape, name="Synthetic C", verbose=True)
    simImageNoised, sigma = createImage(0.3, oObject, psf_A, psf_C, verbose=True)
    print(f"sigma = {sigma}")
    dspNoise = sigma**2 * simImageNoised.shape[0]**2
    dspObjectForEstimation = dspObject[:1024,:1024]
    jointEstimation(oObject, dspObject, psf_A, psf_C, simImageNoised, sigma, xAxis, yAxis, N, "synthetic")
    alpha = marginalEstimation(oObject, dspObject, psf_A, psf_C, simImageNoised, sigma, xAxis, yAxis, N, "synthetic")
    o = estimatedObject(simImageNoised, oObject, dspObjectForEstimation, dspNoise, alpha, psf_A, psf_C, "synthetic")
    dif = oObject.real - o
    saveImplot(dif, "Difference between estimated and real object (synthetic)", filename="Difference between estimated and real object (synthetic)", save=True, plot=True)
    # Qualité d'image comme somme quadratique terme à terme (dB)
    
    
    # Real object / image estimation, with two PSFs
    im, imShape = openImage("extrait1_8b_CD1_20110502165033.tif") # (2*imShape[0],2*imShape[1])
    psf_A, psf_C = updatePsfSize(psf_1, imShape, name="Real A", verbose=True), updatePsfSize(psf_3, imShape, name="Real C", verbose=True)
    # dspObject_est_real = sp.fft.fft2(im) ** 2 # We suppose the object dsp conserved in the image
    # dspObject_est_real = dspObject_est_real[:int(imShape[0]/2),:int(imShape[1]/2)]
    
    rObject, dsprObject = createObject("Real", nPoints=2*imShape[1], sampling_rate=1, p=2.2, k=1, rho_0=0.01)
    dsprObject = dsprObject[:int(imShape[0]),:int(imShape[1])]
    sigmaBis = 0.06 * np.max(im)
    print(f"sigmaBis = {sigmaBis}")
    dspNoise = sigmaBis**2 * imShape[0] * imShape[1]
    # We suppose that the mean value of the object is conserved in the image
    jointEstimation(im, dsprObject, psf_A, psf_C, im, sigmaBis, xAxis, yAxis, N, "real")
    alpha_real = marginalEstimation(im, dsprObject, psf_A, psf_C, im, sigmaBis, xAxis, yAxis, N, "real")
    o = estimatedObject(im, im, dsprObject, dspNoise, alpha_real, psf_A, psf_C, "real")
    
    
    # Synthetic object / image estimation, with three PSFs
    
    # Real object / image estimation, with three PSFs

    
if __name__ == "__main__":
    main()
