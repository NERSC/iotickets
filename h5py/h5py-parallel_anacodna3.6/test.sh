#!/bin/bash
#SBATCH -N 1 #Use 1 node
#SBATCH -t 00:05:00 #Set 5 minute time limit
#SBATCH -q debug #Submit to debug QOS
#SBATCH -C haswell #Use Haswell nodes
module load python/3.6-anaconda-4.4
module load h5py-parallel
srun -n 4 python test.py
