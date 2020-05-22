# TODO:
# Replace suffix rules by pattern rules.
# Include make.sys as dependency.

.SUFFIX.: 
.SUFFIXES.: *.f90

include ./make.sys.gfortran 

modules = mod_util.o mod_types.o mod_MPI_parallel.o
ssetfi_objects =  class_Stack.o \
	          triangular_lattice.o \
	          linked_list_triangular_plaquette.o \
	          SSE_configuration.o \
		  diagonal_update_plaquette.o \
		  quantum_cluster_update_plaquette.o \
		  measurements.o
	  
	
objects = ${modules} ${ssetfi_objects}	  

default: all

all: ssetfi 


ssetfi: ${objects} ssetfi_main.f90
	${F90} ${objects} -o ssetfi ${LFLAGS}

.f90.o:
	${F90} ${FFLAGS} $*.f90

clean:
	rm -f *.mod 
	rm -f *.o
	
clean-dat:
	rm -f *.dat 
	rm -f *.out
	
clean-all: clean clean-dat	
