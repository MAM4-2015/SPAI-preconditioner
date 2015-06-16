# SPAI-preconditioner


_Langage utilisé : Julia_

Les fichiers de ce dépot permettent la réalisation d'un préconditionnement d'une matrice sparse. Très largement inspiré de la méthode de M. J. Grote et T. Huckle _Parallel precontionning with sparse approximate inverses_. 

Fichiers : 

+ Laplacien.jl : 3 fonctions créant la matrice du laplacien 1D, 2D et 3D de taille n fixé par l'utilisateur
+ SPAI.jl : Plusieurs fonctions nécessaires à la construction de l'algorithme
+ SPAIfun.jl : Contient la fonction SPAI qui renvoit la matrice M de préconditionnement, uniquement sur le Laplacien 1D, 2D ou 3D pour d'autre matrice changer l'intitialisation de la matrice A dans la fonction SPAI

Difficulté : 

+ Parallélisation du problème
+ Obligation d'utiliser la fonction _inv_ dans la mise à jour de _mkchap_ car le message stack overflow est renvoyé avec l'utilisation de l'opérateur \

Auteur :

Yacine BERKANE, Kévin ELIE-DIT-COSAQUE, Wahid MAINASSARA, Laurent PAGLIARI

Pour toute question n'hésitez pas à nous contacter.
