import h5py

swmr_file = h5py.File('swmr.h5', 'w', libver='latest')

#This gives an IOError
swmr_file.swmr_mode = True
