CFLAGS= -L/opt/cray/pe/hdf5-parallel/1.8.16/INTEL/15.0/lib
LIBS=-lhdf5_fortran -lhdf5 -lz 
OUTPUTBIN=io_bench
OUTPUTSRC=io_benchmark_blacs.f90
FTN=ftn

.PHONY:all
all:$(OUTPUTBIN)
all:$(OUTPUTBIN)
$(OUTPUTBIN): $(OUTPUTSRC)
	$(FTN) $(OUTPUTSRC)  -o $(OUTPUTBIN) $(CFLAGS) $(LIBS)
.PHONY:clean
clean:
	rm $(OUTPUTBIN)

