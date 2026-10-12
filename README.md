# gf2m

Arithmétique des corps finis GF(2^m) et algèbre linéaire associée, en Rust.

Bibliothèque écrite dans le cadre d'un travail sur la cryptographie post-quantique
fondée sur les codes correcteurs d'erreurs. L'objectif est une implémentation de
référence, destinée à servir à des versions optimisées ultérieures.

## Contenu

- `gf` : éléments de GF(2^m) : addition, multiplication, carré, inverse.
  Le degré `M` et le polynôme irréductible `POLY` sont des paramètres de
  compilation. ( Le degré est limité à M < 64 pour le moment)
- `vect` : vecteurs : somme, produit par un scalaire, produit scalaire.
- `matrice` : somme, produit, transposée, produit matrice-vecteur,
  élimination de Gauss-Jordan, rang, inverse, base du noyau.
