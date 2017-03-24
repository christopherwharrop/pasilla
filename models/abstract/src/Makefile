include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(UTILITIES)/include

all: module_model.o module_config.o

module_model.o: module_model.f90
	$(FC) $(FFLAGS) -c module_model.f90 $(INCLUDES)

module_config.o: module_config.f90
	$(FC) $(FFLAGS) -c module_config.f90 $(INCLUDES)

clean:
	rm -f *.o *.mod

