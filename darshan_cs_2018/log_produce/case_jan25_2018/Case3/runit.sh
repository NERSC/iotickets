#!/bin/bash -l
##SBATCH -p regular
#SBATCH -t 00:20:00
#SBATCH -N 342 #171 #342 #128 #86
#SBATCH -o brtnfld.o%j
##SBATCH -C haswell   #Use Haswell nodes
#Edison has 24 cores per node
cd $SCRATCH
WRKDIR=$SCRATCH/$SLURM_JOB_ID
mkdir $WRKDIR
lfs setstripe -c 64 -S 16m $WRKDIR
cd $WRKDIR
cmd="case1"
tsk=8192 #4096 # -N times 24
cleanup=yes
cp $SLURM_SUBMIT_DIR/$cmd .
srun -n $tsk ./$cmd
cd $SLURM_SUBMIT_DIR
if [ -n "$cleanup" ]; then
  rm -r -f $WRKDIR
fi
