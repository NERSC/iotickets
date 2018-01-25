import h5py
import time
import sys
if len(sys.argv) > 1:
  fname=sys.argv[1]
else:
  fname= 'antielectron_207.2d.h5'
print fname
tstart=time.time()
try:
 fx=h5py.File(fname,'r')
 dx=fx['images']
 tsread=time.time()
 print "starting read slices"
 for i in range(13):
  ttread=time.time()
  temp=dx[i*128:(i+1)*128,0:2,0:240,0:4096]
  ttfinish=time.time()
  print "batch %d costs %f secs"%(i,ttfinish-ttread)
 tstop=time.time()
 print "total read costs %f seconds"%(tstop-tsread)
except Exception as e:
 print "error:%s"%e
 pass
