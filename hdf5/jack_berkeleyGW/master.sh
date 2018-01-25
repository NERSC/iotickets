#!/bin/bash
#SBATCH -N 4
#SBATCH -t 00:30:00
#SBATCH -J blacs
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -C haswell  
#SBATCH -p debug
module load cray-hdf5-parallel/1.8.16
export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_STATS=1
export DARSHAN_DISABLE_SHARED_REDUCTION=1
 
cur=/global/cscratch1/sd/jialin/io-ticket/nerscio/jack_berkeleyGW
#for i in 1 2 4 8 16 32 64 128 244 
for i in 8  
do
  for j in 2 4
  do
   dir=${cur}/ost${i}/
   rm fort*
   rm -rf $dir
   mkdir -p $dir
   echo $dir
   lfs setstripe -c $i -S 1m $dir 
   #cp $cur/io_bench_blacs.x $dir
   npe=32
   srun -N $j -n $((j*npe)) ./io_bench 50000 4 $dir
   ls -alh $dir/*.h5
  done
done
