# TrouNoir

## Compilation
La compilation néscéssite l'installation préalable de BLAS et LAPACK


make compile uniquement la simulation principale

make test compile les tests (placé dans le fichier bin)

make clean supprime les fichier objets et autre fichier intermédiaire de compilation

make mrproper supprime tout les produits de compilation y compris les éxécutables

make debug recompile intégralement avec les options de débuggage

Attention du fait de certaine option d'optimisation les éxécutables ne sont pas portable et le code doit être compilé sur chaque ordinateur sur lesquels il est destiné à être utilisé

## Utilisation
./trouNoir lance le programme principal

./bin/testScurveMesh lance le programme générant les mesh de courbe en S

./bin/testingdicho lance les programmes en rapport avec la dichotomie


physicalsetting contient les paramètres de la simulation

physicalsettingSCurve contient les paramètres pour générer les mesh de courbes en S à l'aide de testScurveMesh
