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
#import line_profiler 
#from memory_profiler import profile 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def create_df(columns, grp, myrank, ranks):
    ''' create a dataframe, given file, group and list of datasets
    arguments:
    columns: datasets to be read as columns in a dataframe
    grp: group name to be read
    fname: file name of the file to read
    myrank: process's rank
    ranks: total number of MPI ranks in the job
    output: a data frame per rank consisting of columns
    '''
    return pd.DataFrame(read_table.read_table(columns, grp, myrank, ranks))

def rename_column(colname, grpname): 
    ''' specific to the cms converted data, need to remove the 
    group name and a period from the data set names so that the 
    column names corresponding to the data sets are shorter and 
    query/filter code is readable. 

    input: data set name and group
    output: a string by removing the group name from the data set names
    '''
    if colname.startswith(grpname.lstrip("/")): 
       colname=colname[len(grpname)+1:]
    return colname

#@profile
def create_df_all(in_dir, grp_name, myrank, ranks): 
    h5files=glob.glob(in_dir+"/*.%s"% 'h5')
    df=pd.DataFrame()
    dfs=[]
    for fname in h5files:
        with h5py.File(fname, 'r', driver='mpio', comm=MPI.COMM_WORLD) as f:
            grp=f["/"+grp_name]
            if (grp_name == 'Info'):
                ds_names=['Info.lumisec', 'Info.runNum', 'Info.evtNum', 'rhoIso']
            else:
                ds_names=grp.keys()
                ds_names.remove('snum')
                ds_names.remove('stype')
            tmp_df = create_df(ds_names, grp, myrank, ranks)
            dfs.append(tmp_df)
    df = pd.concat(dfs)
    del tmp_df
    del dfs
    old_cols=df.columns
    df.columns=[rename_column(e, grp_name) for e in old_cols]
    return df

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: ", sys.argv[0], "input-dir group-name"
        sys.exit(-1)
    s_time = MPI.Wtime()
    in_dir = sys.argv[1]
    grp_name = sys.argv[2]
    myrank = MPI.COMM_WORLD.rank  
    ranks = MPI.COMM_WORLD.size
    elec_df=create_df_all(in_dir, grp_name, myrank, ranks)
    #print len(elec_df)
    info_df=create_df_all(in_dir, "Info", myrank, ranks)
    r_time = MPI.Wtime()
    if myrank==0: print r_time - s_time
    df=pd.merge(elec_df, info_df, how='left', on=['evtNum', 'lumisec', 'runNum'])
    #hist, bin_edges = np.histogram(df['pt'], density=True)
    #hist, bin_edges = np.histogram(df['eta'], density=True)
    #print bin_edges 
    fdf=fe.filter_electrons(df)
    # for pt [  10.00000668   60.96817007  111.93633347  162.90449686  213.87266026
  #264.84082365  315.80898705  366.77715044  417.74531384  468.71347723
  #519.68164062]
    localhistpt, bin_edges = np.histogram(fdf['pt'], bins=[0, 20, 30, 40, 50, 100, 500]) #, density=True) 
    #localhisteta, bin_edges = np.histogram(fdf['eta'], density=True) 
    #print bin_edges
    #local_cnt=len(fdf)
    #global_cnt = MPI.COMM_WORLD.reduce(local_cnt, op=MPI.SUM, root=0)
    globalhist = MPI.COMM_WORLD.reduce(localhistpt, op=MPI.SUM, root=0)
    #print globalhist
    if myrank==0: print globalhist, bin_edges, MPI.Wtime()-r_time
    if myrank==0: 
        #plt.figure()
        plt.hist(globalhist)
        plt.savefig("myplot%s"%ranks) 
       #fig = plt.gcf()
        #plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')
