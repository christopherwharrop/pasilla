GPTL=/contrib/gptl/gptl-v5.4.4_impi_noomp

FC=mpiifort
FCOPT=-O3
FCDEBUG=-g -traceback
FFLAGS=$(FCOPT) -finstrument-functions -mkl=sequential $(FCDEBUG)
INCLUDES=-I$(NETCDF)/include -I$(GPTL)/include
LIBS=-L/contrib/gptl/gptl-v5.4.4_impi_noomp/lib -lgptl

varsolver.exe: varsolver_12aug.o
	$(FC) -o varsolver.exe varsolver_12aug.o $(FFLAGS) $(LIBS)

varsolver_12aug.o: varsolver_12aug.f90
	$(FC) -c varsolver_12aug.f90 $(FFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o

