<!-- <h1 align="center">Wind Watch</h1> -->
<div align="center">
  <img src="/public/Style/wind-watch.jpg">
</div>

<h4 align="center">Web app for computing optimized wind-based sailing courses.</h4>

<!-- <h1 align="center"> </h1> -->

<!-- ![version](https://img.shields.io/badge/version-0.0.1-blueviolet) -->
<div style="flex" align="center">
  <img src="https://img.shields.io/badge/version-1.0.0-blueviolet" alt="version badge">
  <img src="https://img.shields.io/badge/testing-in%20progress-orange" alt="Testing Status">
  <img src="https://img.shields.io/badge/development-in%20progress-orange" alt="Development Status">
  <img src="https://img.shields.io/badge/maintained-yes-brightgreen.svg" alt="Maintenance Status">
  <img src="https://img.shields.io/badge/launched-no-red.svg" alt="Launch Status">
  <img src="https://img.shields.io/badge/license-MIT-blue" alt="License">
</div>

<!-- ![development](https://img.shields.io/badge/development-in%20progress-orange)
![maintenance](https://img.shields.io/badge/maintained-yes-brightgreen.svg)
![launched](https://img.shields.io/badge/launched-no-red.svg)
![License](https://img.shields.io/badge/license-MIT-blue) -->

<br>

---

Myopic deconvolution of retina images acquired by adaptative optics through bayesian estimation.

## Introduction

When an optical system is designed, a theoretical Point Spread Function (PSF) is associated with it. The system built according to this design is calibrated, and the experimental PSF slightly differs from the theoretical one (usually slightly degraded due to the imperfections of the manufacturing processes). For classic usecases of optical systems, like classic photography the in-use PSF will be the experimental one with a good approximation (modulo the aging of the system). However for specific usescases involving middles with variating characteristics (such as optical indices), the PSF is modified with each specific middle. One way of correcting this change of the PSF, due to the changes in the measuring environment, it to introduce adaptative optics (AO) in the design of the optical system. AO will characterize the perturbation (air temperature, winds for astronomical observation involving light propagation in the atmosphere, or of the eye lens and eyeball liquid for retina scanning) and reconstruct the wavefront so that the usecase PSF is partly corrected from this.

This AO design produces better resolved images that can be then improved with digital processing, usually involving deconvolution of the acquired image by the usecase PSF of the system. For deconvolution of images acquired in classic (static) middles, the experimental PSF is used, however here, the estimate of the PSF is to be found to perform such filtering. This situation is called blind deconvolution (or myopic deconvolution since the image still contain information about the PSF). For example, the scheme used in [[1]](#1) uses some available a priori information on the PSF, namely, its positivity and estimates of its ensemble mean and PSD.

In retina imagery, the acquired images are 3 dimensional, however the acquired images are only 2 dimensional (in imagery in general, this concept still holds). So, to look at a specific plane within the depth of focus of the imaging system (outside, the object only contributes by adding background photon flux to the relevant image), it is possible to use the associated PSF. In first approximation, this PSF can be linearly computed between two PSF, in focus and in depth of focus for example. [[2]](#2)

### Goals

- 

### Project structure & description

```shell
.
├── .github                         # GitHub folder
├── .gitignore                      
├── .venv                           # Python virtual environment
├── .vscode                         # VSCode configuration
├── data                            
├── requirements.txt                # Python requeriments
├── README.md                       # README file
│   ├── README.md                   # README describing the data
│   ├── ___.tif                     # Data (images) to analyse
├── outputs                         # Output images
├── Scripts
│   ├── README.md                   # README to get started
│   ├── main.py                     # Next JS App (App Router)
```

## Getting Started

### Requirements

#### Python

- Python 3.12+
- Required Libraries: see `requirements.txt`

### Installation

Run the following command on your local environment:

```shell
python -m venv
source .venv/bin/activate # on unix
.venv/Script/activate # on windows
pip install -r requirements.txt
```

Then run the main function :

```bash
python .Scripts/main.py
```

## Debugging & Testing

Unitarian test are yet to be defined.

## Learn More

Please see wordone.html to learn about the work done in this study and the quotation used from the 2 references [[1]](#1) [[2]](#2).

## Contributions

Everyone is welcome to contribute to this project. Feel free to open an issue if you have any questions or find a bug. Totally open to suggestions and improvements.

### Future Implementation

- Bathymetric maps integration
- Tidal heightsd and current forecast integration
- Boats specification integration
- Loggin methods
- Shipyard implementation : register boats to be used
- Trips list : register trips onto the app for future analytics or sharing with friends

## References

<a id="1">[1]</a>

Jean-Marc Conan, Laurent M. Mugnier, Thierry Fusco, Vincent Michau, and Gerard Rousse | Myopic deconvolution of adaptive optics images by use of object and point-spread function power spectra | APPLIED OPTICS / Vol. 37, No. 21 / 20 July 1998 | © 1998 Optical Society of America

<a id="2">[2]</a>

L. Blanco1, L. M. Mugnier1 |Marginal blind deconvolution of adaptive optics retinal images |7 November 2011 / Vol. 19, No. 23 / OPTICS EXPRESS 23227 |© 2011 OSA

<br>

## License

Licensed under the MIT License, Copyright © 2024

See [LICENSE](LICENSE) for more information.