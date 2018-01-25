#!/bin/bash
#   #SBATCH -p debug
#SBATCH -p regular
#   #SBATCH -d afterok:6121646
#SBATCH -N 30
#SBATCH -C knl,quad,cache
#SBATCH -t 10:00:00
#SBATCH -J bandsq
#   #SBATCH -L project
#SBATCH -S 8
#DW jobdw capacity=5000GB access_mode=striped type=scratch pool=wlm_pool

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#module load espresso/5.4.0

cat run*burst*02*.*
echo "Burst Buffer directory ="$DW_JOB_STRIPED
export TMP_DIR=$DW_JOB_STRIPED
export ESPRESSO_TMPDIR=$DW_JOB_STRIPED

originalDir=$(pwd)
cp -ar pw_in* $DW_JOB_STRIPED 
cp -ar GaN.save $DW_JOB_STRIPED
cp *.UPF $DW_JOB_STRIPED

module load espresso/6.1
srun -n 1800 -c 4 --cpu_bind=cores pw.x -in pw_in.burstBuffer -npool 180 -ntg 1

cp -ar $DW_JOB_STRIPED/GaN.save $originalDir
cp -arf slurm*.out $DW_JOB_STRIPED
du -shc $originalDir/GaN.save


