from collections import defaultdict
from mpi4py import MPI
import glob
import h5py
import math
import numpy as np
import os
import sys

def read_table(columns, grp, myrank, ranks):
    '''Return a dictionary of np.ndarrays read from an HDF group.

    Arguments:
    columns -- a list of the names of columns to be read
    grp -- the name of the group to be read
    fname -- the name of the HDF5 file to be read
    myrank -- the MPI rank of this process
    ranks -- the total number of ranks in the program
    '''
    return {colname: read_ds(colname, grp, myrank, ranks) for colname in columns}

def read_ds(ds_name, grp, myrank, nranks):
    '''Return an np.ndarray read from an HDF dataset.

    Arguments:
    ds -- name of the dataset to be read
    grp -- name of the group in which the dataset is found
    fname -- name of the HDF5 file to be read
    myrank -- the MPI rank of this process
    ranks -- the total number of ranks in the program
    '''
    d = grp[ds_name]
    var_size=d.size
    start = var_size/nranks * myrank
    end = start + int(math.ceil(float(var_size)/nranks))
    if (var_size % nranks) != 0: end = end - 1
    if myrank==nranks-1 : end=var_size
    return d[start:end]

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: ", sys.argv[0], "input-dir group-name"
        sys.exit(-1)

    in_dir = sys.argv[1]
    grp_name = sys.argv[2]
    myrank = MPI.COMM_WORLD.rank  
    nranks = MPI.COMM_WORLD.size
    h5files=glob.glob(in_dir+"/*.%s"% 'hdf5')    
    #h5files=os.listdir(in_dir)
    tables=[]
    d = defaultdict(list)
    s_time = MPI.Wtime()
    for fname in h5files:
        print fname
        with h5py.File(fname, 'r', driver='mpio', comm=MPI.COMM_WORLD) as f:
            grp=f["/"+grp_name]
            ds_names=grp.keys()
            tables.append(read_table(ds_names, grp, myrank, nranks))
    for t in tables:
        for key, value in t.iteritems():
            d[key].append(value)

