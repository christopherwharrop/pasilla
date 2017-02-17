include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(GPTL)/include -I$(UTILITIES)/include $(NETCDF_INCLUDES) -I../../abstract/include -I.
LIBS=-L$(GPTL)/lib -lgptl -L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS)

FFLAGS+=$(GPTLFLAGS)
FFLAGS+=$(ENDIAN_FLAG)

OBJS=module_QG_model.o nag.o

all: QG.exe

QG.exe: runqgmodel.f90 $(OBJS)
	$(FC) $(FFLAGS) -o QG.exe runqgmodel.f90 $(OBJS) $(INCLUDES) $(LIBS)

nag.o: nag.f
	$(FC) $(FFLAGS) -c nag.f

module_QG_model.o: module_QG_model.f90 nag.o
	$(FC) $(FFLAGS) -c module_QG_model.f90 $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod *.a

