include ../../macros.mk

INCLUDES=-I$(GPTL)/include -I../../utilities/include -I../../models/abstract/include -I../../models/sine/include -I.
LIBS=-L$(GPTL)/lib -lgptl -L../../models/sine/lib -lsine -L$(NETCDF)/lib -lnetcdf -lnetcdff
FFLAGS+=-openmp

varsolver.exe: varsolver.f90 module_varsolver.o
	$(FC) -o varsolver.exe varsolver.f90 module_varsolver.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver.o: module_varsolver.f90
	$(FC) -c module_varsolver.f90 $(FFLAGS) $(GPTLFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod

