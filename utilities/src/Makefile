include ../../macros.mk

OBJS=module_kind.o module_constants.o

libpasilla_util.a: $(OBJS)
	ar rcs libpasilla_util.a $(OBJS)

module_kind.o: module_kind.f90
	$(FC) $(FFLAGS) -c module_kind.f90

module_constants.o: module_constants.f90
	$(FC) $(FFLAGS) -c module_constants.f90

clean:
	rm -f *.o *.a *.mod

