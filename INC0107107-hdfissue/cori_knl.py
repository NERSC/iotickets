#!/usr/bin/env python3.4
import os
from math import *



# -- Function to create submission script --------------------------------------

def gen_script(filename):

   # On Cori KNL, nodes have 96GB of memory and 68 physical cores. We subtract 4GB
   # for the OS.  Also NERSC suggests setting aside 2 or 4 cores (we use 4) for the
   # OS, so there are 92 GB and 64 cores available.
   mem_per_core = 92/64.  # GB

   # Cores are grouped into "tiles" of 2 cores with shared L2 cache, so we prefer
   # an even number of cores per task.
   if mem_per_task <= mem_per_core:
      phys_cores_per_task = 1
   else:
      phys_cores_per_task = int(ceil(mem_per_task/(2.0*mem_per_core))) * 2

   # There are 4 hyperthreads (logical cores) per physical core.
   logc_cores_per_task = phys_cores_per_task*4

   # Testing indicates the hyperthreading hurts overall, so we just use 1 thread
   # per physical core.
   omp_num_threads = phys_cores_per_task

   tasks_per_node = int(floor(64./phys_cores_per_task))
   nodes = int(ceil(total_tasks/float(tasks_per_node)))

   omp_num_threads=min(omp_num_threads,logc_cores_per_task)

   s="""#!/bin/bash
#SBATCH -N {nodes}
#SBATCH -C knl,quad,cache
#SBATCH -p {queue}
#SBATCH -J {name}
#SBATCH -S 4
#SBATCH --mail-user=jalnliu@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH --qos=premium

# OpenMP settings:
export OMP_NUM_THREADS={omp_num_threads}
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Job info
echo "input file: {infile}"
echo "number of MPI tasks: {total_tasks}"
echo "mem per task: {mem_per_task} GB"
echo "physical cores per task: {phys_cores_per_task}"
echo "logical cores per task: {logc_cores_per_task}"
echo "tasks per node: {tasks_per_node}"
echo "number of nodes: {nodes}"
echo "OpenMP threads per task: {omp_num_threads}"
echo "OMP_PLACES: $OMP_PLACES"


# Check process and thread affinity
srun -n {total_tasks} -c {logc_cores_per_task} --cpu_bind=cores check-hybrid.intel.cori

# Run the application
srun -n {total_tasks} -c {logc_cores_per_task} --cpu_bind=cores {exe} mc {infile}
"""

   # Use global and local variables to format the script
   d = locals().copy(); d.update(globals())
   script = s.format_map(d)

   print("Creating " + filename)
   myfile = open(filename,"w")
   myfile.write(script)
   myfile.close()



# -- Function to create input file ---------------------------------------------


def gen_input(filename):
   s="""# Output files
result file: {name}.res
data file: {name}.h5

# System
lattice points: 13
g: -44.044
nspecies: 2
spin: F
spherical cutoff: F

n1: 65
n2: 65
min n1:64
max n1:66
min n2:64
max n2:66

dbeta: {dbeta}
beta: {beta}
SO: 0.0

# Observables
ODLRO: T
Nk: T
heat capacity: T
dT/T: 0.02

# Sampling
samples: 25600
sweeps per sample: 5
thermalization samples: 0
dsigma: 0.6
random seed: {seed}
"""

   d = locals().copy(); d.update(globals())
   s1 = s.format_map(d)

   print("Creating " + filename)
   myfile = open(filename,"w")
   myfile.write(s1)
   myfile.close()



# -- Generate inputs and scripts -----------------------------------------------

seed=21231
for db in [128]:
   for beta in [0.9375, 1.00, 1.078125, 1.140625]:
#   for beta in [1.078125, 1.140625]:
      total_tasks =  50      # Number of MPI tasks (processes)
      mem_per_task = 5      # GB  (ADJUST ME)

      seed=seed+5
      dbeta = 1./db
      name = "lx13.nspc2.SO0.db{db}.b{beta}".format(db=db,beta=beta)
      exe="lafmc"
      queue="regular"            # "debug" or "regular"

      infile = name + ".in"
      gen_input(infile)

      scriptname = name + ".sh"
      gen_script(scriptname)

      # Submit
      os.system("sbatch " + scriptname)
