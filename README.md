# SPAI-preconditioner


_Langage utilisé : Julia_

Les fichiers de ce dépot permettent la réalisation d'un préconditionnement d'une matrice sparse. Très largement inspiré de la méthode de M. J. Grote et T. Huckle _Parallel preconditionning with sparse approximate inverses_. 

Fichiers : 

+ Laplacien.jl : 3 fonctions créant la matrice du laplacien 1D, 2D et 3D de taille n fixé par l'utilisateur
+ SPAI.jl : Plusieurs fonctions nécessaires à la construction de l'algorithme
+ SPAIfun.jl : Contient la fonction SPAI qui renvoit la matrice M de préconditionnement, uniquement sur le Laplacien 1D, 2D ou 3D pour d'autre matrice changer l'intitialisation de la matrice A dans la fonction SPAI
+ _SPAIfun2.jl : Contient la fonction SPAI2. Il s'agit d'une copie de la fonction SPAI de SPAIfun.jl mais avec rajout de certaines conditions pour éviter des erreurs sur des cas très particuliers. La fonction SPAI devrait être suffisante pour tous les types de matrices mais par soucis de précaution nous préférons mettre cette deuxième fonction pour être sur de pouvoir gérer tous les cas_

Difficulté : 

+ Parallélisation du problème
+ Obligation d'utiliser la fonction _inv_ dans la mise à jour de _mkchap_ car le message stack overflow est renvoyé avec l'utilisation de l'opérateur \

Auteur :

Yacine BERKANE, Kévin ELIE-DIT-COSAQUE, Wahid MAINASSARA, Laurent PAGLIARI

Pour toute question n'hésitez pas à nous contacter.

Professeur encadrant :

Thierry CLOPEAU
