#!/bin/bash -l
#SBATCH --account=m1523
#SBATCH -p debug
#SBATCH -N 4
#SBATCH -t 00:07:00
#SBATCH -J my_job
#SBATCH -o test_128_4_2.o
#SBATCH -L SCRATCH
#SBATCH -C haswell  #Use Haswell nodes

module unload python
#module load python
module load python/2.7-anaconda
#module load h5py-parallel/2.6.0
#module load mpi4py
module load h5py-parallel/2.7.1
#module load mpi4py
module list

#export OMP_NUM_THREADS=4
#srun -n 2 -c 4 python demo1.py
#srun -n 64 python ~/spark-hdf5-cms/python/create_df.py
#srun -n 128 python create_df.py /global/cscratch1/sd/ssehrish/h5out_merged/ Electron 
srun -n 128 python-mpi hist_test.py /global/cscratch1/sd/ssehrish/h5out_merged/ Electron 0 haswell
#srun -n 128 python ~/hep_hpc_python/read_table.py /global/cscratch1/sd/ssehrish/lariat-hdf5-1/Run3a/ rawDigits_u  
