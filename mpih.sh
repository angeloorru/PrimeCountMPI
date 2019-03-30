#!/bin/bash
#PBS -N my_mpi_job
#PBS -l nodes=1:ppn=4

cd $PBS_O_WORKDIR
mpirun -np 4 -machinefile $PBS_NODEFILE ./primeCount >& ./primeCount.log

