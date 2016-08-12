#FC=ifort
FC=ifort
FFLAGS=-O3 -g -traceback
INCLUDES=-I. -I$(NETCDF)/include
LIBS=-L$(NETCDF)/lib -lnetcdf -lnetcdff

lorenz.exe: module_lorenz96.o lorenz96.f90
	$(FC) $(FFLAGS) -o lorenz.exe lorenz96.f90 module_lorenz96.o $(LIBS)

module_kind.o: module_kind.f90
	$(FC) $(FFLAGS) -c module_kind.f90

module_lorenz96.o: module_lorenz96.f90 module_kind.o
	$(FC) $(FFLAGS) -c module_lorenz96.f90 $(INCLUDES)

clean:
	rm -f lorenz.exe *.o *.mod

