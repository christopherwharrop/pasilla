include ../../macros.mk

INCLUDES=-I$(GPTL)/include -I.
LIBS=-L$(GPTL)/lib -lgptl

varsolver.exe: varsolver.f90 module_varsolver.o
	$(FC) -o varsolver.exe varsolver.f90 module_varsolver.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver.o: module_varsolver.f90
	$(FC) -c module_varsolver.f90 $(FFLAGS) $(GPTLFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod

