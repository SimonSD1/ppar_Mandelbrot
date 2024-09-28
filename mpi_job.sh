#!/bin/bash

# Fichier où stocker les sorties
output_file="mpi_output.txt"

# Nettoyer le fichier de sortie au début
> $output_file

# Boucle sur les différentes valeurs de n (nombre de processus)
for n in $(seq 1 1 6); do
    echo "Running mpiexec with --n $n" >> $output_file
    # Exécuter mpiexec avec la valeur de n et rediriger la sortie vers le fichier
    mpiexec --n $n  ./mandel 12000 12000 0.35 0.355 0.353 0.358 1000 >> $output_file 2>&1
    echo "Finished mpiexec with --n $n" >> $output_file
done
