#!/bin/bash
#SBATCH -p debug
#SBATCH -N 3
#SBATCH -t 00:30:00
#SBATCH -J Sp-S1-64
#SBATCH -e mysparkjob_64-core_%j-debug.err
#SBATCH -o mysparkjob_64-core_%j-debug.out
#SBATCH --ccm
#SBATCH --qos=premium
#module unload spark/hist-server
module load spark
#module load collectl
#start-collectl.sh 
start-all.sh

export LD_LIBRARY_PATH=$LD_LBRARY_PATH:$PWD/lib

###load single large hdf5 file####
repartition="1000"
inputfile="/global/cscratch1/sd/dbin/ost-248/udf/fake-2d-tiny.h5p"
app_name="SparkS1"
dataset="/testg/testd"

spark-submit --verbose\
  --master $SPARKURL\
  --name $app_name \
  --driver-memory 100G\
  --executor-cores 32 \
  --driver-cores 32  \
  --num-executors=2 \
  --executor-memory 105G\
  --class org.nersc.io.readtest\
  --conf spark.eventLog.enabled=true\
  --conf spark.eventLog.dir=$SCRATCH/spark/spark_event_logs\
  target/scala-2.10/h5spark-assembly-1.0.jar \
  $repartition $inputfile $dataset 


##rm /global/cscratch1/sd/jialin/spark_tmp_dir/*
stop-all.sh
#stop-collectl.sh
