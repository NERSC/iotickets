#!/bin/bash
#SBATCH -N 4
#SBATCH -C knl,quad,cache
#SBATCH -p regular
#SBATCH -J lx13.nspc2.SO0.db128.b1.140625
#SBATCH -S 4
#SBATCH --mail-user=jalnliu@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH --qos=premium

# OpenMP settings:
export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Job info
echo "input file: lx13.nspc2.SO0.db128.b1.140625.in"
echo "number of MPI tasks: 50"
echo "mem per task: 5 GB"
echo "physical cores per task: 4"
echo "logical cores per task: 16"
echo "tasks per node: 16"
echo "number of nodes: 4"
echo "OpenMP threads per task: 4"
echo "OMP_PLACES: $OMP_PLACES"


# Check process and thread affinity
srun -n 50 -c 16 --cpu_bind=cores check-hybrid.intel.cori

# Run the application
srun -n 50 -c 16 --cpu_bind=cores lafmc mc lx13.nspc2.SO0.db128.b1.140625.in
