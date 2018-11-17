#!/bin/bash
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:35:00
#SBATCH -A dasrepo
#DW jobdw capacity=10GB access_mode=striped type=scratch

mkdir $DW_JOB_STRIPED/inputdir
cp $SCRATCH/test.h5 $DW_JOB_STRIPED/inputdir
file1=$DW_JOB_STRIPED/inputdir/test.h5
file2=/project/projectdirs/das/jialin/test_copy.h5
srun -n 1 python test.py $file1 $file2 

