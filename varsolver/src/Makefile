include ../../macros.mk

INCLUDES=-I$(GPTL)/include -I../../utilities/include -I../../models/abstract/include -I../../models/sine/include -I../../models/lorenz96/include -I.
LIBS=-L$(GPTL)/lib -lgptl -L../../models/sine/lib -lsine -L$(NETCDF)/lib -lnetcdf -lnetcdff -L../../models/lorenz96/lib -llorenz96
FFLAGS+=$(OPENMP)
FFLAGS+=$(GPTLFLAGS)


all: varsolver_sine.exe varsolver_l96.exe

varsolver_sine.exe: varsolver_sine.f90 module_varsolver_sine.o
	$(FC) -o varsolver_sine.exe varsolver_sine.f90 module_varsolver_sine.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver_sine.o: module_varsolver_sine.f90
	$(FC) -c module_varsolver_sine.f90 $(FFLAGS) $(INCLUDES)

varsolver_l96.exe: varsolver_l96.f90 module_varsolver_l96.o
	$(FC) -o varsolver_l96.exe varsolver_l96.f90 module_varsolver_l96.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver_l96.o: module_varsolver_l96.f90
	$(FC) -c module_varsolver_l96.f90 $(FFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod

