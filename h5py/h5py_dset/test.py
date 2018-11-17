import h5py
import time
fx=h5py.File('haha_nofill.h5','w') 
spaceid = h5py.h5s.create_simple((30000,8000000)) 
plist = h5py.h5p.create(h5py.h5p.DATASET_CREATE) 
plist.set_fill_time(h5py.h5d.FILL_TIME_NEVER) 
plist.set_chunk((60, 15625)) 
start=time.time()
datasetid =h5py.h5d.create(fx.id,"data",h5py.h5t.NATIVE_FLOAT, spaceid, plist) 
dset = h5py.Dataset(datasetid) # then you can use normal h5py api to deal with this dset object 
fx.close() 
end=time.time()
print ('%f seconds'%(end-start))
