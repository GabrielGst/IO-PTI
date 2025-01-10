# IO-PTI
Myopic deconvolution of retina images acquired by adaptative optics through bayesian estimation.

## Context
When an optical system is designed, a theoretical Point Spread Function (PSF) is associated with it. The system built according to this design is calibrated, and the experimental PSF slightly differs from the theoretical one (usually slightly degraded due to the imperfections of the manufacturing processes). For classic usecases of optical systems, like classic photography the in-use PSF will be the experimental one with a good approximation (modulo the aging of the system). However for specific usescases involving middles with variating characteristics (such as optical indices), the PSF is modified with each specific middle. One way of correcting this change of the PSF, due to the changes in the measuring environment, it to introduce adaptative optics (AO) in the design of the optical system. AO will characterize the perturbation (air temperature, winds for astronomical observation involving light propagation in the atmosphere, or of the eye lens and eyeball liquid for retina scanning) and reconstruct the wavefront so that the usecase PSF is partly corrected from this.

This AO design produces better resolved images that can be then improved with digital processing, usually involving deconvolution of the acquired image by the usecase PSF of the system. For deconvolution of images acquired in classic (static) middles, the experimental PSF is used, however here, the estimate of the PSF is to be found to perform such filtering. This situation is called blind deconvolution (or myopic deconvolution since the image still contain information about the PSF). For example, the scheme used in [[1]](#1) uses some available a priori information on the PSF, namely, its positivity and estimates of its ensemble mean and PSD.

In retina imagery, the acquired images are 3 dimensional, however the acquired images are only 2 dimensional (in imagery in general, this concept still holds). So, to look at a specific plane within the depth of focus of the imaging system (outside, the object only contributes by adding background photon flux to the relevant image), it is possible to use the associated PSF. In first approximation, this PSF can be linearly computed between two PSF, in focus and in depth of focus for example.

## Bayesian estimation
In astronomy, the estimated objects are the object and the PSF, while in retina imagery, these are the object and the parameter $\alpha$ which in the ends defines the PSF through linear fit. In astronomy, joint estimation is used in conjunction with 

### Joint estimation

### Marginal estimation


## References
<a id="1">[1]</a> 
Jean-Marc Conan, Laurent M. Mugnier, Thierry Fusco, Vincent Michau, and Gerard Rousse | 
Myopic deconvolution of adaptive optics images by use of object and point-spread function power spectra | 
APPLIED OPTICS / Vol. 37, No. 21 / 20 July 1998 | 
© 1998 Optical Society of America

