GPTL=/contrib/gptl/gptl-v5.5_nompi_noomp

FC=ifort
FCOPT=-O3
FCDEBUG=-g -traceback
GPTLFLAGS=-finstrument-functions 
FFLAGS=$(FCOPT) -mkl=parallel $(FCDEBUG) -openmp 

INCLUDES=-I$(GPTL)/include -I.
LIBS=-L$(GPTL)/lib -lgptl

varsolver.exe: varsolver.f90 module_varsolver.o
	$(FC) -o varsolver.exe varsolver.f90 module_varsolver.o $(FFLAGS) $(INCLUDES) $(LIBS)

module_varsolver.o: module_varsolver.f90
	$(FC) -c module_varsolver.f90 $(GPTLFLAGS) $(FFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod

