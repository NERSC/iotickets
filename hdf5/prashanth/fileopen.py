# Test to multiply two HDF5 datasets

from mpi4py import MPI
import numpy as np
import h5py
import time
import sys

# filename_x = '/global/project/projectdirs/dasrepo/BoLBO/x_reduced_2.hdf5'
filename_x = '/global/project/projectdirs/dasrepo/BoLBO/x_test.hdf5'
f_x = h5py.File(filename_x,'r')

# filename = '/global/project/projectdirs/dasrepo/BoLBO/x_transpose_test.hdf5'
#filename = '/global/cscratch1/sd/pselv/x_transpose/x_transpose_test.hdf5'
filename = '/global/cscratch1/sd/jialin/hdf-data/prashan/x_transpose_test.hdf5'
colw = 0

#length_x = f_x['x'].shape[1]
length_x = 2*50
length_y = f_x['x'].shape[0]

comm = MPI.COMM_WORLD
nproc = comm.Get_size()

h5f = h5py.File(filename, 'w', driver='mpio', comm=MPI.COMM_WORLD)
rank = comm.Get_rank()

dset = h5f.create_dataset('x_trans', (length_x,length_y), dtype='f32')
h5f.atomic = False
total_write=length_x*length_y*32/1024/1024/1024 #GB
each_write=total_write/nproc*1024 #MB

length_rank=length_x / nproc
length_last_rank=length_x -length_rank*(nproc-1)
comm.Barrier()
timestart=MPI.Wtime()
start=int(rank*length_rank)
end=start+int(length_rank)
if rank==nproc-1: #last rank
    end=start+int(length_last_rank)

print "Process number: %d, Start index: %d, Finish index: %d" %(rank,start,end)
if rank==0:
        print "start transposing"
temp = f_x['x'][:,start:end]
#temp =  np.transpose(f_x['x'][:,start:end])
timetrans=MPI.Wtime()
if rank==0:
         print "finish transposing"
         print "Transpose time %f" %(timetrans-timestart)

if rank==0:
	print "start writing"
#if colw==1:
#	with dset.collective:
#		dset[start:end,:] = temp
#else:
#	dset[start:end,:] = temp
comm.Barrier()

timeend=MPI.Wtime()
if rank==0:
    print "finish writing"
    if colw==1:
    	print "Collective write time %f" %(timeend-timestart)
    else :
        #print "Transpose time %f" %(timetrans-timestart)
	print "Independent write time %f" %(timeend-timetrans)
    print "Data size x: %d y: %d" %(length_x, length_y)
    print "File size ~%d GB" % (length_x*length_y/1024.0/1024.0/1024.0*4)
    print "Number of processes %d" %nproc
h5f.close()
f_x.close()
