GPTL=/contrib/gptl/gptl-v5.4.4_impi_noomp

FC=mpiifort
FCOPT=-O3
FCDEBUG=-g -traceback
FFLAGS=$(FCOPT) -finstrument-functions -mkl=sequential $(FCDEBUG)
INCLUDES=-I$(NETCDF)/include -I$(GPTL)/include
LIBS=-L/contrib/gptl/gptl-v5.4.4_impi_noomp/lib -lgptl

varsolver.exe: varsolver.o
	$(FC) -o varsolver.exe varsolver.o $(FFLAGS) $(LIBS)

varsolver.o: varsolver.f90
	$(FC) -c varsolver.f90 $(FFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o

