include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(UTILITIES)/include $(NETCDF_INCLUDES) -I../../abstract/include
LIBS=-L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS)

FFLAGS+=$(OPENMP)

all: sine.exe sineTL_test.exe sineADJ_test.exe libsine.a

libsine.a: $(OBJS)
	ar rcs libsine.a module_sine.o

sine.exe: module_sine.o sine.f90
	$(FC) $(FFLAGS) -o sine.exe sine.f90 module_sine.o ../../abstract/src/module_model.o $(INCLUDES) $(LIBS)

sineTL_test.exe: module_sine.o sineTL_test.f90
	$(FC) $(FFLAGS) -o sineTL_test.exe sineTL_test.f90 module_sine.o ../../abstract/src/module_model.o $(INCLUDES) $(LIBS)

sineADJ_test.exe: module_sine.o sineADJ_test.f90
	$(FC) $(FFLAGS) -o sineADJ_test.exe sineADJ_test.f90 module_sine.o ../../abstract/src/module_model.o $(INCLUDES) $(LIBS)

module_sine.o: module_sine.f90
	$(FC) $(FFLAGS) -c module_sine.f90 $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod *.a

