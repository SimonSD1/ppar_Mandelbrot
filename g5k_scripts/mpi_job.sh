#!/bin/bash

output_file="output_temps_dans_le_c.txt"

>$output_file

echo "Contenu de \$OAR_NODEFILE:" >>$output_file
cat $OAR_NODEFILE >>$output_file
echo "" >>$output_file 

function run_mpi() {
    local n=$1
    local total_time=0

    echo "Running mpiexec with --n $n" >>$output_file

    for i in $(seq 1 5); do

        mpiexec --n $n --mca pml ^ucx --hostfile $OAR_NODEFILE  ./mandel 12000 12000 0.35 0.355 0.353 0.358 1000 >> $output_file 2>&1

    done
}

for n in $(seq 1 1 10); do
    run_mpi $n
done

for n in $(seq 20 10 150); do
    run_mpi $n
done
