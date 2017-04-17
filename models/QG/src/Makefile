include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(GPTL)/include -I$(UTILITIES)/include $(NETCDF_INCLUDES) -I../../abstract/include -I.
LIBS=-L$(GPTL)/lib -lgptl -L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS)

FFLAGS+=$(GPTLFLAGS)
FFLAGS+=$(ENDIAN_FLAG)

OBJS=module_QG_model.o module_QG_config.o module_QG_writer.o module_QG_reader.o module_QG_ggsp.o nag.o

all: QG.exe libqg.a

libqg.a: $(OBJS) # module_qg_tl.o module_qg_adj.o
	ar rcs libqg.a $(OBJS) # module_qg_tl.o module_qg_adj.o

QG.exe: QG.f90 $(OBJS)
	$(FC) $(FFLAGS) -o QG.exe QG.f90 $(OBJS) $(INCLUDES) $(LIBS)

nag.o: nag.f
	$(FC) $(FFLAGS) -c nag.f

module_QG_model.o: module_QG_model.f90 module_QG_config.o module_QG_ggsp.o nag.o
	$(FC) $(FFLAGS) -c module_QG_model.f90 $(INCLUDES)

module_QG_config.o: module_QG_config.f90
	$(FC) $(FFLAGS) -c module_QG_config.f90 $(INCLUDES)

module_QG_writer.o: module_QG_writer.f90 module_QG_model.o
	$(FC) $(FFLAGS) -c module_QG_writer.f90 $(INCLUDES)

module_QG_reader.o: module_QG_reader.f90 module_QG_model.o
	$(FC) $(FFLAGS) -c module_QG_reader.f90 $(INCLUDES)

module_QG_ggsp.o: module_QG_ggsp.f90 nag.o
	$(FC) $(FFLAGS) -c module_QG_ggsp.f90 $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod *.a

