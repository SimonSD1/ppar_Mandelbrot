# Mandelbrot

Premier TP de parallélisation du calcul de l'ensemble de Mandelbrot et analyse des performances.

Compilation : mpicc mandel.c rasterfile.h -o mandel -lm -O3

Exemple d'éxecution : mpiexec ./mandel 12000 12000 0.35 0.355 0.353 0.358 1000
