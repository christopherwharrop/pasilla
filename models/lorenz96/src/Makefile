include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(UTILITIES)/include -I$(NETCDF)/include
LIBS=-L$(UTILITIES)/lib -lpasilla_util -L$(NETCDF)/lib -lnetcdf -lnetcdff

lorenz96.exe: module_lorenz96.o lorenz96.f90
	$(FC) $(FFLAGS) -o lorenz96.exe lorenz96.f90 module_lorenz96.o $(INCLUDES) $(LIBS)

module_lorenz96.o: module_lorenz96.f90
	$(FC) $(FFLAGS) -c module_lorenz96.f90 $(INCLUDES)

clean:
	rm -f lorenz96.exe *.o *.mod

