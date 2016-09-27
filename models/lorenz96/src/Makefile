include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(UTILITIES)/include $(NETCDF_INCLUDES)
LIBS=-L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS)

all: lorenz96.exe lorenz96TL_test.exe lorenz96ADJ_test.exe

lorenz96.exe: module_lorenz96.o lorenz96.f90
	$(FC) $(FFLAGS) -o lorenz96.exe lorenz96.f90 module_lorenz96.o $(INCLUDES) $(LIBS)

lorenz96TL_test.exe: module_lorenz96.o lorenz96TL_test.f90
	$(FC) $(FFLAGS) -o lorenz96TL_test.exe lorenz96TL_test.f90 module_lorenz96.o $(INCLUDES) $(LIBS)

lorenz96ADJ_test.exe: module_lorenz96.o lorenz96ADJ_test.f90
	$(FC) $(FFLAGS) -o lorenz96ADJ_test.exe lorenz96ADJ_test.f90 module_lorenz96.o $(INCLUDES) $(LIBS)

module_lorenz96.o: module_lorenz96.f90
	$(FC) $(FFLAGS) -c module_lorenz96.f90 $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod

