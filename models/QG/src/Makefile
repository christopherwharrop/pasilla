include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(UTILITIES)/include $(NETCDF_INCLUDES) -I../../abstract/include -I.
LIBS=-L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS)

FFLAGS+=$(GPTLFLAGS)
FFLAGS+=$(ENDIAN_FLAG)

OBJS=comqg.o qgmodel.o nag.o

all: QG.exe

QG.exe: runqgmodel.f90 $(OBJS)
	$(FC) $(FFLAGS) -o QG.exe runqgmodel.f90 $(OBJS) $(INCLUDES) $(LIBS)

nag.o: nag.f
	$(FC) $(FFLAGS) -c nag.f

comqg.o: comqg.f90
	$(FC) $(FFLAGS) -c comqg.f90

qgmodel.o: qgmodel.f90 nag.o comqg.o
	$(FC) $(FFLAGS) -c qgmodel.f90

clean:
	rm -f *.exe *.o *.mod *.a

