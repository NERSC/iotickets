from mpi4py import MPI
import glob
import h5py
import numpy as np
import os
import pandas as pd
import sys
import read_table
import tables
import filter_electrons as fe
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import create_df as cd

if __name__ == "__main__":
  if len(sys.argv) < 5:
    print "Usage: ", sys.argv[0], "input-dir group-name run-num platform"
    sys.exit(-1)
  myrank = MPI.COMM_WORLD.rank  
  ranks = MPI.COMM_WORLD.size
  #if myrank==0: print "Run", "Nodes", "Cores", "Label", "Time"
  for x in range(0, 1):
    s_time = MPI.Wtime()
    in_dir = sys.argv[1]
    grp_name = sys.argv[2]
    run_num = sys.argv[3]
    platform = sys.argv[4]
    elec_df = cd.create_df_all(in_dir, grp_name, myrank, ranks)
    info_df = cd.create_df_all(in_dir, "Info", myrank, ranks)
    r_time = MPI.Wtime()
    if myrank==0: print ranks/32, ranks, "Read", r_time - s_time, run_num, platform 
    s_time = MPI.Wtime()
    df = pd.merge(elec_df, info_df, how='left', on=['evtNum', 'lumisec', 'runNum'])
    if myrank==0: print ranks/32, ranks, "Merge", MPI.Wtime() - s_time, run_num, platform
    s_time = MPI.Wtime()
    fdf = fe.filter_electrons(df)
    if myrank==0: print ranks/32, ranks, "Filter", MPI.Wtime() - s_time, run_num, platform
    s_time = MPI.Wtime()
    localhistpt, bin_edges = np.histogram(fdf['pt'], bins=[0, 20, 30, 40, 50, 100, 500]) 
    if myrank==0: print ranks/32, ranks, "Myplot", MPI.Wtime() - s_time, run_num, platform
    s_time = MPI.Wtime()
    globalhist = MPI.COMM_WORLD.reduce(localhistpt, op=MPI.SUM, root=0)
    if myrank==0: print ranks/32, ranks, "Plot", MPI.Wtime() - s_time, run_num, platform
    if myrank==0: 
        plt.hist(globalhist)
        plt.savefig("myplot%s"%ranks)
