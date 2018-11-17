from mpi4py import MPI
import h5py
mpi_rank = MPI.COMM_WORLD.Get_rank()
#import h5py

fx=h5py.File('output.h5','w',driver='mpio',comm=MPI.COMM_WORLD)
print(mpi_rank)
fx.close()
