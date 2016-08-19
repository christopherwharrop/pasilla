UTILITIES=../../utilities

FC=ifort
FCFLAGS=-O3 -g -traceback
INCLUDES=-I$(UTILITIES)/include -I$(NETCDF)/include
LIBS=-L$(UTILITIES)/lib -lpasilla_util -L$(NETCDF)/lib -lnetcdf -lnetcdff

lorenz96.exe: module_lorenz96.o lorenz96.f90
	$(FC) $(FCFLAGS) -o lorenz96.exe lorenz96.f90 module_lorenz96.o $(INCLUDES) $(LIBS)

module_lorenz96.o: module_lorenz96.f90
	$(FC) $(FCFLAGS) -c module_lorenz96.f90 $(INCLUDES)

clean:
	rm -f lorenz96.exe *.o *.mod

