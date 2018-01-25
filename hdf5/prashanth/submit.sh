#!/bin/sh -l

#SBATCH -N 1
#SBATCH -p debug
#SBATCH -t 00:12:00
#SBATCH -J transpose
#SBATCH -e %j.err
#SBATCH -o %j.out
#SBATCH -V

cd $SLURM_SUBMIT_DIR
module unload python
module unload altd
module swap PrgEnv-intel PrgEnv-gnu
module load python
module load mpi4py
module load h5py-parallel
srun -n 2 python-mpi fileopen.py
