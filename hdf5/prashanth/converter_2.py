#Row Matrix Transposing to Column Matrix
#Input is a list of 2D netcdf file
#Output is one single HDF5 transposed file

from mpi4py import MPI
from netCDF4 import Dataset
import h5py
import numpy as np
from os import listdir
from os.path import isfile, join
import sys
import time
import math, sys

rank = MPI.COMM_WORLD.Get_rank()
mpi_info = MPI.Info.Create()
numProcs = MPI.COMM_WORLD.Get_size()

numProcessesPerNode = 10

procslist = np.arange(numProcs)

def status(message, ranks=procslist):
    if rank in ranks:
        print "%s, process %d: %s" % (time.asctime(time.localtime()), rank, message)

def report(message):
    status(message, [0])

def reportbarrier(message):
    MPI.COMM_WORLD.Barrier()
    report(message)

# maps the chunkidx (0...numlevdivs=numWriters) to the rank of the process that should write it out
def chunkidxToWriter(chunkidx):
   return (chunkidx*numProcessesPerNode)%numProcs # should always write from different nodes
datapath = "/global/cscratch1/sd/gittens/large-climate-dataset/data"
filelist = [fname for fname in listdir(datapath) if fname.endswith(".nc")]
#filelist = [fname for fname in filelist[:200]]
varname = (sys.argv[1])
report("Using %d processes" % numProcs)
report("Writing variable %s" % varname)
report("Found %d input files, starting to open" % len(filelist))

assert( len(filelist) % numProcs == 0)

# open all the files associated with this process and keep handles around (avoid metadata costs)
myfiles = [fname for (index, fname) in enumerate(filelist) if (index % numProcs == rank)]
numFilesPerProc = len(myfiles)
myhandles = [None]*len(myfiles)
for (idx, fname) in enumerate(myfiles):
    myhandles[idx] = Dataset(join(datapath, fname), "r") 

reportbarrier("Finished opening all files")

#varnames = ["T", "U", "V", "Q", "Z3"]
#varnames = ["T"]
varnames = [varname]
numvars = len(varnames)
numtimeslices = 8
numlevels = 30
numlats = 768
numlongs = 1152
numlevdivs = 256
flattenedlength = numlevels*numlats*numlongs
numRows = flattenedlength*numvars
numCols = len(filelist)*numtimeslices
rowChunkSize = numlats*numlongs/numlevdivs

numWriters = numlevdivs 
coloutname = "production/colfilenames" + varname 
foutname = "production/" + varname
assert ((numlats * numlongs) % numlevdivs == 0)

report("Creating files and datasets")

# Write out the names of the files corresponding to each column in the final output matrix, ordered 0,...,numCols
# simplify!
colfilenames = MPI.COMM_WORLD.gather(myfiles, root=0)
if rank == 0:
   colnames = [item for sublist in colfilenames for item in sublist] 
   with open(coloutname, "w") as colfout:
       np.save(colfout, np.array(colnames))

propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
propfaid.set_fapl_mpio(MPI.COMM_WORLD, mpi_info)
fid = h5py.h5f.create(join(datapath, foutname + ".h5"), flags=h5py.h5f.ACC_TRUNC, fapl=propfaid)
fout = h5py.File(fid)

# Don't use filling 
spaceid = h5py.h5s.create_simple((numvars*numlevels*numlats*numlongs, numCols))
plist = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
plist.set_fill_time(h5py.h5d.FILL_TIME_NEVER)
datasetid = h5py.h5d.create(fout.id, "rows", h5py.h5t.NATIVE_FLOAT, spaceid, plist)
rows = h5py.Dataset(datasetid)
reportbarrier("Finished creating files and datasets")

localcolumncount = numFilesPerProc*numtimeslices
curlevdata = np.empty((numlats*numlongs, localcolumncount), \
        dtype=np.float32)
chunktotransfer = np.empty((rowChunkSize*localcolumncount,), dtype=np.float32)

listwriter = map(chunkidxToWriter, np.arange(numWriters))
if rank in listwriter:
    collectedchunk = np.ascontiguousarray(np.empty((numCols*rowChunkSize,), \
            dtype=np.float32))
    chunktowrite = np.ascontiguousarray(np.empty((rowChunkSize, numCols), \
            dtype=np.float32))
else:
    collectedchunk = None
curlevdatatemp=np.ascontiguousarray(np.zeros((numlats*numlongs*numtimeslices), \
            dtype=np.float32))

for (varidx,curvar) in enumerate(varnames): 
    reportbarrier("Writing variable %d/%d: %s" % (varidx + 1, numvars, curvar))

    for curlev in np.arange(numlevels):

        # load the data for this level from my files
        reportbarrier("Loading data for level %d/%d" % (curlev + 1, numlevels))
        for (fhidx, fh) in enumerate(myhandles):
            if fh[curvar].shape[0] < numtimeslices and fh[curvar].shape[0] >0:
                status("File %s has only %d timesteps for variable %s, simply repeating the first timestep" % (myfiles[fhidx], fh[curvar].shape[0], curvar))
                for idx in np.arange(numtimeslices):
                    curlevdatatemp[numlats*numlongs*idx:numlats*numlongs*(idx+1)] = fh[curvar][0, curlev, ...].flatten()
		curlevdata[:, fhidx*numtimeslices: (fhidx + 1)*numtimeslices]=curlevdatatemp.reshape(numlats*numlongs, numtimeslices)
 		curlevdatatemp=np.ascontiguousarray(np.zeros((numlats*numlongs*numtimeslices), dtype=np.float32)) 
            elif fh[curvar].shape[0] ==0:
		status("File %s has only %d timesteps for variable %s, simply repeating the first timestep" % (myfiles[fhidx], fh[curvar].shape[0], curvar))
		curlevdata[:, fhidx*numtimeslices: (fhidx + 1)*numtimeslices] = \
			curlevdatatemp.reshape(numlats*numlongs, numtimeslices)
	    else:
                curlevdata[:, fhidx*numtimeslices: (fhidx + 1)*numtimeslices] = \
                    fh[curvar][:, curlev, ...].reshape(numlats*numlongs, numtimeslices)
	
        reportbarrier("Done loading data for this level")
        
        # write out this level in several chunks of rows
        reportbarrier("Gathering data for this level from processes to writers")
        for chunkidx in np.arange(numlevdivs):
            startrow = chunkidx*rowChunkSize
            endrow = startrow + rowChunkSize
            chunktotransfer[:] = curlevdata[startrow:endrow, :].flatten()
            MPI.COMM_WORLD.Gather(chunktotransfer, collectedchunk, root = chunkidxToWriter(chunkidx))
	reportbarrier("Done gathering")
	
	reportbarrier("Writing data for this level on writers")
        for chunkidx in np.arange(numlevdivs):
            if rank == chunkidxToWriter(chunkidx):
                for processnum in np.arange(numProcs):
                    startcol = processnum*localcolumncount
                    endcol = (processnum+1)*localcolumncount
                    startidx = processnum*(localcolumncount * rowChunkSize)
                    endidx = (processnum + 1)*(localcolumncount *rowChunkSize)
                    chunktowrite[:, startcol:endcol] = np.reshape(collectedchunk[startidx:endidx], \
                            (rowChunkSize, localcolumncount))	
		start_rows = varidx*numlevels*numlats*numlongs+curlev*numlats*numlongs+rowChunkSize*chunkidx
		end_rows = start_rows+rowChunkSize
		rows[start_rows:end_rows, :] = chunktowrite 
        reportbarrier("Done writing")

for fh in myhandles:
    fh.close()
fout.close()

