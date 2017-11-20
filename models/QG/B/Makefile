include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(GPTL)/include -I$(UTILITIES)/include $(NETCDF_INCLUDES) -I../../abstract/include -I.
LIBS=-L$(GPTL)/lib -lgptl -L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS)

all: bmatrix.exe latlon.exe sqrt_b.exe

%.exe:%.f90
	$(FC) -o $@ $(FFLAGS) $(INCLUDES) $(LIBS) $<

clean:
	rm -f *.exe

