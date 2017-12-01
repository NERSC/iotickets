import h5py
from mpi4py import MPI
import numpy as np
import time
import sys
comm = MPI.COMM_WORLD
nproc = comm.Get_size()
fx=h5py.File('antielectron_207.2d.h5','r',driver='mpio',comm=MPI.COMM_WORLD)
rank=comm.Get_rank()
dx=fx['images']
comm.Barrier()
timestart=MPI.Wtime()
temp=np.empty(dx.shape,dx.dtype)
d_shape=dx.shape
d_type=dx.dtype
temp[rank*128:(rank+1)*128,:]=dx[rank*128:(rank+1)*128,0:2,0:240,0:4096]
comm.Barrier()
timeend=MPI.Wtime()
if rank==0:
  print "number of ranks %d"%nproc
  print "total read cost %f"%(timeend-timestart)
#print "rank %d start %d, value %f"%(rank,rank*128+1,dx[rank*128+1,0,0,0])
fx.close()
#start writing out the raw data
comm.Barrier()
timestart=MPI.Wtime()
fx_dec=h5py.File('antielectron_207.2d_raw.h5','a',driver='mpio',comm=MPI.COMM_WORLD)
dx_dec=fx_dec.create_dataset('images',d_shape,dtype='f2')
dx_dec[rank*128:(rank+1)*128,:]=temp[rank*128:(rank+1)*128]
fx_dec.close()
comm.Barrier()
timeend=MPI.Wtime()
if rank==0:
  print "total write cost %f"%(timeend-timestart)
