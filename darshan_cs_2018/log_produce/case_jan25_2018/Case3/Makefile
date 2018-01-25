HOST := $(shell hostname)

ifeq ($(findstring edison, $(HOST)),edison)
 CC = cc
 FC = ftn
else
 CC = mpicc
 FC = mpif90
endif

F90FLAGS = -O3

OBJF90 = case1.o

all:  case1

case1: $(OBJ) $(OBJF90)
	$(FC) $(F90FLAGS) $(OBJF90) -o $@

.SUFFIXES: .o .F90

.F90.o:
	$(FC) $(F90FLAGS) -c $< -o $@

clean:
	rm -f  ../NonContig/*.o *.o *.mod case1 timing datafile.mpio

spotless:
	rm -f ../NonContig/*.o *.o *.mod case1 *~ timing datafile.mpio
