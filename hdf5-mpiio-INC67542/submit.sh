#!/bin/bash
#SBATCH -p debug 
#SBATCH -N 2
#SBATCH -t 00:10:00
#SBATCH -J hdf5writeTest_2
#SBATCH -e %j.err
#SBATCH -o %j.out


srun -n 64 ./hdf5writeTest_2 -x 64
