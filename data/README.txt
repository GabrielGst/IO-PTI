Le fichier fits contient 3 PSFs, qui peuvent �tre utilis�es pour simuler des images :

- PSF_1 = d�foc nulle,

- PSF_2 = d�foc RMS  pi/2

- PSF_3 = defoc RMS  pi.

Commencer par n'utiliser que les 2 PSFs extremes, pour pouvoir *tracer* le crit�re en fonction d'un alpha scalaire sans utiliser d'optimiseur :

psf = alpha*psf_1 + (1-alpha)*psf_3, avec par exemple alpha=0.2.


Utiliser un objet conforme � l'a priori ie une realisation d'un champ gaussien avec DSP donn�e par le mod�le utilis� :

pente. avec p=2 ca converge + lentement

p = 1 ; 
rho = toutes les coordonn�es des fr�quences spatiales en polaire

modele de DSP objet parametrique
dspo = k/((rho/rho0)^p+1) ; 

generation du champ aleatoire
tfo = randomn(seed, NP, NP) * sqrt(dspo)

TF inverse du champ
objet = partie r�elle de la fft inverse de tfo monodimensionnelle
objet = objet-min(objet)

afficher l'objet