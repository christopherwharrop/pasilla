include ../../macros.mk

OBJS=module_kind.o module_constants.o module_netcdf_utilities.o

libpasilla_util.a: $(OBJS)
	ar rcs libpasilla_util.a $(OBJS)

module_kind.o: module_kind.f90
	$(FC) $(FFLAGS) -c module_kind.f90

module_constants.o: module_constants.f90 module_kind.o
	$(FC) $(FFLAGS) -c module_constants.f90

module_netcdf_utilities.o: module_netcdf_utilities.f90 module_kind.o
	$(FC) $(FFLAGS) -c module_netcdf_utilities.f90 $(NETCDF_INCLUDES)

clean:
	rm -f *.o *.a *.mod

