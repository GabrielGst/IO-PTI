Le fichier fits contient 3 PSFs, qui peuvent être utilisées pour simuler des images :

- PSF_1 = défoc nulle,

- PSF_2 = défoc RMS  pi/2

- PSF_3 = defoc RMS  pi.

Commencer par n'utiliser que les 2 PSFs extremes, pour pouvoir *tracer* le critère en fonction d'un alpha scalaire sans utiliser d'optimiseur :

psf = alpha*psf_1 + (1-alpha)*psf_3, avec par exemple alpha=0.2.


Utiliser un objet conforme à l'a priori ie une realisation d'un champ gaussien avec DSP donnée par le modèle utilisé :

pente. avec p=2 ca converge + lentement

p = 1 ; 
rho = toutes les coordonnées des fréquences spatiales en polaire

modele de DSP objet parametrique
dspo = k/((rho/rho0)^p+1) ; 

generation du champ aleatoire
tfo = randomn(seed, NP, NP) * sqrt(dspo)

TF inverse du champ
objet = partie réelle de la fft inverse de tfo monodimensionnelle
objet = objet-min(objet)

afficher l'objet