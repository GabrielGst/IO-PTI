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

Here, we want to estimate the object using the inverse convolution of the image and the Point Spread Function (PSF): as an image can be described as the convolution of the object and the PSF, understood as the Fourier Transform of the wave incident on the pupil plane, it is only logical to compute the object that way.

$$
i = o * PSF
$$

where 

$$
h = PSF = \mid FT[\Phi(u) \, P(u)] \mid ^2
$$

We will call Optical Transfer Function (OTF) the Fourier transform of the PSF :

$$
\tilde{h} = OTF = FT[PSF]
$$

However, the PSF changes according to the distance of the optical system to the object. We thus also have to estimate the right PSF using two extrema : the in-focus and out-of-focus PSF. We thus define the right PSF as :

$$
PSF = \alpha \times PSF_{In} + (1-\alpha) \times PSF_{Out}
$$

In reality, the image contains noise that needs to be accounted for in the formula, leading to further consideration when computing inverse convolution and impacting the expression of such filters.

### Retina imaging

$$ \mathbf{i}_{3 \mathrm{D}}=\mathbf{h}_{3 \mathrm{D}} *_{3 \mathrm{D}} \mathbf{o}_{3 \mathrm{D}}+\mathbf{n} $$

$$ \mathbf{i}_{2 \mathrm{D}}=\mathbf{h}_{2 \mathrm{D}} *_{2 \mathrm{D}} \mathbf{o}_{2 \mathrm{D}}+\mathbf{n} $$

We assume that our object is shift invariant along the optical axis:
$$
o_{3 \mathrm{D}}(x, y, z)=o_{2 \mathrm{D}}(x, y) \alpha(z)
$$
where $\alpha(z)$ is the normalized flux emitted by the plane at depth $z$ (considering $\int \alpha(z) \mathrm{d} z=1$).

$$
h_{2 \mathrm{D}}(x, y) \approx \sum_j \alpha_j h_j(x, y)
$$
with $h_j(x, y) \triangleq h_{3 \mathrm{D}}\left(x, y, z_j\right)$ the 2D lateral PSF at depth $z_j$ and $\alpha_j=\alpha\left(z_j\right) \Delta z_j$ where $\Delta z_j$ is the effective thickness of the $j$ th layer. We define $\alpha=\left\{\alpha_j\right\}_j$ as the vector of unknowns that parameterize the PSF. $\alpha$ is normalized $\left(\Sigma \alpha_j=1\right)$ and each parameter is positive $\left(\alpha_j \geq 0\right)$.

We search for $h_{2 \mathrm{D}}$ as a linear combination of a basis of PSF's, each corresponding to a given plane. In the following, we consider short-exposure diffractive PSF's so that each $h_j$ can be computed from the residual aberrations measured with a WFS and the knowledge of the defocus of plane $z_j$.


Hereafter are extractions from [[1]](#1) and [[2]](#2) that describes bayesian principles in joint and marginal estimation.

### Bayesian estimation [[1]](#1)

In stochastic approaches the object is seen as one realization of a stochastic process. The object is endowed with an a priori distribution $p(\mathbf{o})$, and Bayes' rule combines the likelihood of the data $p(\mathbf{i} \mid \mathbf{o})$ with this a priori distribution into the a posteriori probability distribution $p(\mathbf{o} \mid \mathbf{i})$ :
$$
p(\mathbf{o} \mid \mathbf{i}) \propto p(\mathbf{i} \mid \mathbf{o}) \;p(\mathbf{o})
$$

This leads to two commonly used object estimation methods: the MAP estimation and the MMSE estimation. On the one hand, the MAP estimation defines the restored object as the most probable object, given the data:
$$
\hat{\mathbf{o}}_{\text {map }}=\underset{\mathbf{o}}{\arg \max } \; p(\mathbf{o} \mid \mathbf{i})
$$

On the other hand, the MMSE estimator is defined as the one that minimizes, on average, the distance with the true object:
$$
\hat{\mathbf{o}}_{\mathrm{mmse}}=\underset{\hat{o}}{\arg \min } \; E\left(\|\hat{\mathbf{o}}-\mathbf{o}\|^2\right)
$$
where $E()$ stands for the mathematical expectation with respect to the object and to the image noise. It can be shown that this estimator is the mean object with respect to the a posteriori probability distribution ${ }^{23,24}$ :
$$
\hat{\mathbf{o}}_{\text {mmse }}=E(\mathbf{o} \mid \mathbf{i})=\int \mathbf{o} \, p(\mathbf{o} \mid \mathbf{i}) \, \mathrm{d} \mathbf{o}
$$

In general, the calculation of the MMSE estimator is not tractable unless the estimator is assumed to be linear. This assumption leads to the Wiener filter. It is important to note that in the case of joint Gaussian statistics for the noise and the object, the Wiener, the MMSE, and the MAP estimators are identical. 

<!-- ${ }^{23}$ -->

Choosing the regularization function consist in finding the right model for the regularization parameter, ie the PSD of the *a priori* distribution (of the object here). The maximization criterion is composed of two terms classicaly, the likelihood terms, usually a least square term, and the regularization function.

$$
\begin{aligned} p(\mathbf{o} \mid \mathbf{i}) & \propto p(\mathbf{i} \mid \mathbf{o}) p(\mathbf{o}) \\ & \propto \exp \left[-1 / 2(\mathbf{i}-H \mathbf{o})^t R_n^{-1}(\mathbf{i}-H \mathbf{o})\right] \\ & \times \exp \left[-1 / 2\left(\mathbf{o}-\mathbf{o}_{\mathbf{m}}\right)^t R_o^{-1}\left(\mathbf{o}-\mathbf{o}_{\mathbf{m}}\right)\right],\end{aligned}
$$

### Myopic deconvolution [[2]](#2)

We therefore generalized the deconvolution scheme to the case of myopic deconvolution, in which both the object and the PSF have to be restored. Similarly to what was done for o, the PSF can be considered a stochastic process. Since the PSF can be considered the temporal average of a large number of short-exposure PSF's, its a priori statistics can reasonably be assumed to be Gaussian, and the estimator becomes

$$
\begin{aligned}
{[\hat{\mathbf{o}}, \hat{\mathbf{h}}] } & =\underset{\mathbf{o}, \mathbf{h}}{\arg \max } \; p(\mathbf{o}, \mathbf{h} \mid \mathbf{i}) \\
& =\underset{\mathbf{o}, \mathbf{h}}{\arg \max } \; p(\mathbf{i}|\mathbf{o}, \mathbf{h}) \, p(\mathbf{o}) \, p(\mathbf{h}) \\
& =\underset{\mathbf{o}, \mathbf{h}}{\arg \min } \; J(\mathbf{o}, \mathbf{h})
\end{aligned}
$$

with a new criterion $J(\mathbf{o}, \mathbf{h})$, which is now a function of $\boldsymbol{o}$ and $\mathbf{h}$. This criterion has three terms: one is the opposite of the log likelihood of the data, one is an object regularization term, and one is a PSF regularization term, similar to a recently suggested deterministic approach. For stationary Gaussian noise, this criterion can be written as

$$
J(\mathbf{o}, \mathbf{h}) =  \sum_f\left[\frac{|\tilde{\mathbf{h}}(f) \tilde{\mathbf{o}}(f)-\tilde{\mathbf{i}}(f)|^2}{\operatorname{PSD}_n(f)}+\frac{\left|\tilde{\mathbf{o}}(f)-\tilde{\mathbf{o}}_{\mathbf{m}}(f)\right|^2}{\operatorname{PSD}_o(f)}+\frac{\left|\tilde{\mathbf{h}}(f)-\tilde{\mathbf{h}}_{\mathbf{m}}(f)\right|^2}{\operatorname{PSD}_h(f)}\right]
$$

<!-- $$
\begin{aligned}
J(\mathbf{o}, \mathbf{h})= & \sum_f\left[\frac{|\tilde{\mathbf{h}}(f) \tilde{\mathbf{o}}(f)-\tilde{\mathbf{i}}(f)|^2}{\operatorname{PSD}_n(f)}+\frac{\left|\tilde{\mathbf{o}}(f)-\tilde{\mathbf{o}}_{\mathbf{m}}(f)\right|^2}{\operatorname{PSD}_o(f)}\right. \\
& \left.+\frac{\left|\tilde{\mathbf{h}}(f)-\tilde{\mathbf{h}}_{\mathbf{m}}(f)\right|^2}{\operatorname{PSD}_h(f)}\right]
\end{aligned}
$$ -->

where $\operatorname{PSD}_h$ is the spatial PSD of the PSF, and $\tilde{\mathbf{h}}_{\mathrm{m}}$ is the ensemble mean OTF (Fourier transform of the ensemble mean PSF). Again, when the noise is not Gaussian (which is the case in astronomical imaging), this estimator is not a true MAP estimator but a myopic RLS estimator, unless the first term of the criterion is replaced with the $\log$ probability of the noise.

The last term (regularization on the PSF) cannot be ignored, otherwise the myopic deconvolution usually leads to the trivial solution: a Dirac function for the PSF and an object equal to the image. $\mathrm{PSD}_h$ is expressed simply as a function of the first two moments of the OTF:
$$
\operatorname{PSD}_h(f)=E\left[\left|\tilde{\mathbf{h}}(f)-\tilde{\mathbf{h}}_{\mathrm{m}}(f)\right|^2\right]=E\left[|\tilde{\mathbf{h}}(f)|^2\right]-\left|\tilde{\mathbf{h}}_{\mathrm{m}}(f)\right|^2
$$

The restoration quality can be quantitatively evaluated by the calculation of a distance to the true object $\mathbf{o}$, defined in [[2]](#2) as :
$$ d(\hat{\mathbf{o}}, \mathbf{o})=\left[\frac{1}{N_{\text {pix }}} \sum_{\text {pixels }}|\hat{\mathbf{o}}(r)-\mathbf{o}(r)|^2\right]^{1 / 2}(photons/pixel) $$

In astronomy, the estimated objects are the object and the PSF, while in retina imagery, these are the object and the parameter $\alpha$ which in the ends defines the PSF through linear fit. In astronomy, joint estimation is used in conjunction with a positivity constrain. In both cases, the image is reconstructed and filtered with an *a posteriori* statistics.

Hereafter we switch back to the estimation of $\alpha$ that leads to estimating $h$.

### Joint estimation

$$
\begin{aligned}
{[\hat{\mathbf{o}}, \hat{\mathbf{\alpha}}] } & =\underset{\mathbf{o}, \mathbf{\alpha}}{\arg \max } \; p(\mathbf{o}, \mathbf{\alpha} \mid \mathbf{i} ; \mathbf{\theta}) \\
& =\underset{\mathbf{o}, \mathbf{h}}{\arg \max } \; p(\mathbf{i} \mid \mathbf{o}, \mathbf{\alpha}; \mathbf{\theta}) \, p(\mathbf{o;\theta}) \, p(\mathbf{\alpha;\theta})
\end{aligned}
$$

where $p( \mathbf{o}, \alpha \mid \mathbf{i} ; \boldsymbol{\theta})$ is the joint probability density of the data, of the 2D object, and of the PSF decomposition coefficients ($\alpha$). It may depend on set of regularization parameters or hyperparameters ($\theta$). $p(\mathbf{i} \mid \mathbf{o}, \alpha ; \theta)$ is the likelihood of the data, $p(\mathbf{o} ; \theta)$ is the a priori probability density function of the object $\mathbf{o}$ and $p(\alpha ; \theta)$ is the a priori probability density function of the coefficients $\alpha$. In the following, we will not use any regularization on the set of coefficients $\alpha$ because we do not have any probability law for the PSF coefficients. However, since we only need to estimate a small number of these coefficients, this is not a problem.

The noise on the images is mainly photon noise which has a Poisson distribution. However, AO retinal images are dominated by a strong and quite homogeneous background. In the following, we will therefore assume that the noise is stationary white Gaussian with a variance $\sigma^2$. For the object, we choose a stationary Gaussian prior probability distribution with a mean value $\mathbf{o}_{\mathrm{m}}$ and a covariance matrix $\mathbf{R}_0$. The set of hyperarameters is therefore $\theta=\left(\sigma^2, \mathbf{o}_m, \mathbf{R}_o\right)$. Under these assumptions, we have:
$$
\begin{aligned}
& p(\mathbf{i}, \mathbf{o}, \alpha ; \theta)=\frac{1}{(2 \pi)^{\frac{N^2}{2}} \sigma^{N^2}} \exp \left(-\frac{1}{2 \sigma^2}(\mathbf{i}-\mathbf{H o})^t(\mathbf{i}-\mathbf{H o})\right) \\
& \times \frac{1}{(2 \pi)^{\frac{N^2}{2}} \operatorname{det}\left(\mathbf{R}_{\mathrm{o}}\right)^{1 / 2}} \exp \left(-\frac{1}{2}\left(\mathbf{o}-\mathbf{o}_{\mathrm{m}}\right)^t \mathbf{R}_{\mathrm{o}}^{-1}\left(\mathbf{o}-\mathbf{o}_{\mathrm{m}}\right)\right),
\end{aligned}
$$

$$
\hat{\mathbf{o}}(\alpha, \theta)=\left(\mathbf{H}^t \mathbf{H}+\sigma^2 \mathbf{R}_0^{-1}\right)^{-1}\left(\mathbf{H}^t \mathbf{i}+\sigma^2 \mathbf{R}_0^{-1} \mathbf{o}_{\mathrm{m}}\right)
$$

Since the matrices $\mathbf{H}$ (convolution operator) and $\mathbf{R}_0$ (covariance matrix of an object with a stationary probability density) are Toeplitz-block-Toeplitz, we can write the joint criterion $J_{\text {jmap }}$ and the analytical expression of the object $\hat{\boldsymbol{o}}(\alpha, \theta)$ in the Fourier domain with a circulant approximation:

$$
\hat{\tilde{\mathbf{o}}}(\alpha)=\frac{\tilde{h}^*(v) \tilde{i}(v)+\frac{S_{\mathrm{n}}}{S_0(v)} \tilde{o}_{\mathrm{m}}(v)}{|\tilde{h}(v)|^2+\frac{S_0}{S_0(v)}}
$$

where $S_{\mathrm{n}}$ is the noise power spectral density (PSD), $S_{\mathrm{o}}$ is the object PSD (the new set of hyperparameters in the Fourier domain is $\left\{S_{\mathrm{n}}, S_{\mathrm{o}}\right\}$ ), $v$ is the spatial frequency and $\tilde{x}$ denotes the two-dimensional Fast Fourier Transform of $x$. $\hat{\tilde{\mathbf{o}}}(\alpha)$ is the estimated object after classical Wiener filtering of the image $\mathbf{i}$ and is easily computed. If we substitute $\hat{\tilde{\mathbf{o}}}(\alpha)$ in $J_{\text {jmap }}$, we obtain a new expression $J'_{\text {jmap }}$ that does not depend explicitly on the object:

$$
\begin{aligned}
& J_{\text {jmap }}^{\prime}(\alpha)=\frac{1}{2} N^2 \ln S_{\mathrm{n}}+\frac{1}{2} \sum_v \ln S_{\mathrm{o}}(v) \\
& \quad+\frac{1}{2} \sum_v \frac{1}{S_{\mathrm{o}}(v)} \frac{\left|\tilde{i}(v)-\tilde{h}(v) \tilde{o}_{\mathrm{m}}(v)\right|^2}{|\tilde{h}(v)|^2+\frac{S_0}{S_{\mathrm{o}}(v)}}
\end{aligned}
$$

The joint MAP solution is thus the pair $(\hat{o}(\alpha), \alpha)$ for the value of $\alpha$ that minimizes Eq. 13 .

How can we find $\alpha$ that minimizes $J_{\text {jmap }}$ ?

A simulated image is built in the following manner:
$$
\mathbf{i}=\left(\alpha * \mathbf{h}_{\mathrm{foc}}+(1-\alpha) \mathbf{h}_{\mathrm{defoc}}\right) * \mathbf{o}+\mathbf{n}
$$

The defocus is equal to $\pi$ radian RMS, the noise n is stationary, gaussian with a standard deviation $\sigma = 0.01 \times max(o)$, corresponding roughly to photon noise for an average of 10 000 photons/pixel, and $\alpha =  0.3$.

We assume for the sake of this simulation that the object PSD $S_{\mathrm{o}}$ and the noise PSD $S_{\mathrm{n}}$ are known although it is not the case in practice. Here, we therefore perform a so-called "supervised" estimation of $\alpha$ : we compute the joint criterion $J_{\text {jmap }}\left(\alpha ; S_{\mathrm{o}}, S_{\mathrm{n}}\right)$ for values of $\alpha$ ranging from 0 to 1 to find the value of $\alpha$ that minimizes the joint criterion, and according to an object DSP :
 
$$
DSP_o = S_o = \frac{k}{1+(\rho/\rho_0)^p}
$$

with $(k, p, \rho_0) = (1, 3, 1)$. Please note that in practice, the estimation of the object DSP should be estimated from the acquired image. For example, the expression used here can be used with estimation of the set of parameters $(k, p, \rho_0)$.

### Marginal estimation

The principle of marginal estimation is to integrate the object $o$ out of the problem (i.e., marginalize the posterior likelihood). We integrate the joint probability of the object $o$ and the PSF parameters $\alpha$ over all the possible values of object $o$.
$$
\hat{\alpha}=\underset{\alpha}{\arg \max } \int p( \mathbf{o}, \alpha \mid \mathbf{i} ; \theta) \, \mathrm{d} \mathbf{o}
$$

Marginalization reduces the number of unknowns to be retrieved (from the total number of pixels of the image + the PSF parameters in the joint estimation case to just a few PSF parameters) and gives us a true maximum likelihood or maximum a posteriori (depending on the prior on the estimated parameters) estimator of the parameters of interest (namely, the PSF parameters). After estimation of the PSF parameters $\alpha$, the object is restored by Wiener filtering of the image with the estimated global PSF and hyperparameters.

$$
\hat{\alpha}_{ML}=\underset{\alpha}{\operatorname{argmax}} \; p( \alpha \mid \mathbf{i} ; \theta)=\underset{\alpha}{\operatorname{argmax}} \; p(\mathbf{i} \mid \alpha ; \theta) \, p(\alpha ; \theta)
$$

We keep the assumptions made for the joint estimation: a stationary white Gaussian noise with variance $\sigma^2$, stationary Gaussian prior probability distribution with a mean value $0_{\mathrm{m}}$ and covariance matrix $R_0$ for the object. Since $i$ is a linear combination of a Gaussian object and a Gaussian noise, it is also Gaussian. Its associated probability density reads:
$$
p(\mathbf{i} \mid \alpha ; \theta)=A\left(\operatorname{det} \mathbf{R}_{\mathrm{i}}\right)^{-1 / 2} \exp \left(-\frac{1}{2}\left(\mathbf{i}-\mathbf{i}_{\mathrm{m}}\right)^t \mathbf{R}_{\mathrm{i}}^{-1}\left(\mathbf{i}-\mathbf{i}_{\mathrm{m}}\right)\right)
$$

where $A$ is a constant, $R_i$ is the image covariance matrix and $\mathrm{i}_m=\mathrm{H} \mathrm{o}_{\mathrm{m}}$. Since we only need to estimate a small number of parameters, there is no need to regularize the solution over $\alpha$. We therefore use a Maximum Likelihood (ML) estimator rather than a Maximum A Posteriori (MAP) estimator. Maximizing $p(\mathrm{i} \mid \alpha ; \theta)$ is equivalent to minimizing the opposite of its logarithm:

$$
J_{\mathrm{ML}}(\alpha)=J_{\mathrm{jmap}}^{\prime}(\alpha)+\frac{1}{2} \sum_v \ln \left(|\tilde{h}(v)|^2+\frac{S_{\mathrm{n}}}{S_0(v)}\right)-\frac{1}{2} N^2 \ln S_{\mathrm{n}}
$$

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