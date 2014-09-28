# Default build type

BUILD = debug
#BUILD = production

# Default executable

Exe = comm3d

# Compiler options for different machines

HOST = $(shell hostname)
UNAME = $(shell uname)

# Mac (using gfortran and openmpi from MacPorts)

ifeq ($(findstring Darwin,$(UNAME)),Darwin)
 F95 = openmpif90 -fdefault-real-8 -fdefault-double-8
 LD = openmpif90
 ifeq ($(BUILD),debug)
  F95FLAGS = -g -Wall -Wextra -Wconversion -fbounds-check -fbacktrace \
	-fimplicit-none -std=f2008
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -O5 -Wuninitialized
 endif
 LDFLAGS = $(F95FLAGS)
 LIBS = 
 INCL = 
endif	

# Stanford CEES

ifeq ($(findstring cees,$(HOST)),cees)
 F95 = mpif90 -r8 -i4 -i_dynamic
 LD = mpif90
 ifeq ($(BUILD),debug)
  F95FLAGS =  -g -check all -check noarg_temp_created -warn all \
	-traceback -ftrapuv -fpe0 -fp-stack-check -fltconsistency -std03
 endif
 ifeq ($(BUILD),profile)
  F95FLAGS = -g -pg
 endif
 ifeq ($(BUILD),production)
  F95FLAGS = -g -O2
 endif
 LDFLAGS = $(F95FLAGS)
 LIBS = 
 INCL = 
endif

# Files

Files = comm3d_total.f90 mpi3dbasic.f90 mpi3dcomm.f90 mpi3dio.f90 mpi3d_interface.f90 mpi3d_interface_test.f90 fruit.f90 fruit_util.f90 

Obs = 	$(Files:.f90=.o)

all: $(Exe)
	mv *.o *.mod obj/ 

$(Exe): $(Obs)
	$(LD) $(LDFLAGS) -o $(Exe) \
	$(Obs) $(LIBS)

%.o : %.f90
	$(F95) $(F95FLAGS) -c $< -o $@ $(INCL)


clean:
	rm -f obj/*.o obj/*.mod *.o *.mod *.dat $(Exe)

.SUFFIXES: .o .f90

# DO NOT DELETE THIS LINE - used by make depend
comm3d_total.o: mpi3dbasic.o mpi3dcomm.o mpi3dio.o mpi3d_interface.o mpi3d_interface_test.o

mpi3dcomm.o: mpi3dbasic.o

mpi3dio.o: mpi3dbasic.o

mpi3d_interface_test.o: mpi3d_interface.o fruit.o mpi3dcomm.o 

mpi3d_interface.o: mpi3dcomm.o mpi3dbasic.o


fruit.o: fruit_util.o

