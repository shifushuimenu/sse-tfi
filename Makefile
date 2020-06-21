# TODO:
# Replace suffix rules by pattern rules.
# Include make.sys as dependency for all object files.

.SUFFIXES:            # delete all default suffixes 
.SUFFIXES: .f90 .o    # add suffixes 

include ./make.sys.gfortran 

modules = mod_types.o mod_util.o mod_MPI_parallel.o
ssetfi_objects =  class_Stack.o \
		  SSE_configuration.o \
	      triangular_lattice.o \
	      linked_list_triangular_plaquette.o \
	      diagonal_update_plaquette.o \
	      ssetfi_globals.o \
		  test_helper.o \
		  quantum_cluster_update_plaquette.o \
		  measurements.o
	  
	
objects = ${modules} ${ssetfi_objects}	  

default: all

all: ssetfi 


ssetfi: ${objects} ssetfi_main.f90 make.sys.gfortran
	${F90} ${LFLAGS} ${objects} -o ssetfi ssetfi_main.f90

test_suite: ${objects} test_suite.f90 
	${F90} ${LFLAGS} ${objects} -o test test_suite.f90 
.f90.o: 
	${F90} ${FFLAGS} -c $*.f90

clean:
	rm -f *.mod 
	rm -f *.o
	
clean-dat:
	rm -f *.dat 
	rm -f *.out
	
clean-all: clean clean-dat	
