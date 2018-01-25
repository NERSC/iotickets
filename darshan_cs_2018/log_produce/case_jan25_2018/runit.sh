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
out="2"
cmd="case1"
tsk=8192 #4096 # -N times 24
cleanup=yes
cp $SLURM_SUBMIT_DIR/$cmd .
# find number of MPI tasks per node
#set TPN=`echo $SLURM_TASKS_PER_NODE | cut -f 1 -d \(`
# find number of CPU cores per node
#set PPN=`echo $SLURM_JOB_CPUS_PER_NODE | cut -f 1 -d \(`
#echo " $tsk $TPN $PPN $SLURM_TASKS_PER_NODE $SLURM_JOB_CPUS_PER_NODE"
rm -f $SLURM_SUBMIT_DIR/${out}_$tsk
for i in `seq 1 4`; do
  srun -n $tsk ./$cmd
done
cp timing $SLURM_SUBMIT_DIR/timing_${tsk}_${out}
cd $SLURM_SUBMIT_DIR
if [ -n "$cleanup" ]; then
  rm -r -f $WRKDIR
fi
