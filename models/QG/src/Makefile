include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(GPTL)/include -I$(UTILITIES)/include $(NETCDF_INCLUDES) -I../../abstract/include -I.
LIBS=-L$(GPTL)/lib -lgptl -L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS) -L tapenade/ADFirstAidKit -lstack

FFLAGS+=$(GPTLFLAGS)
FFLAGS+=$(ENDIAN_FLAG)
FFLAGS+=-g -traceback -O0

#OBJS=module_QG_model.o module_QG_config.o module_QG_writer.o module_QG_reader.o module_QG_ggsp.o nag.o nagb.o module_QG_tl.o module_QG_adj.o module_QG_util.o module_QGTL_writer.o module_QGADJ_writer.o
OBJS=module_QG_model.o module_QG_config.o module_QG_writer.o module_QG_reader.o module_QG_ggsp.o nag.o module_QG_tl.o module_QG_adj.o module_QG_util.o module_QGTL_writer.o module_QGADJ_writer.o

all: QG.exe QGTLAD.exe QGTL_test.exe QGADJ_test.exe libqg.a

libqg.a: $(OBJS) module_QG_tl.o module_QG_adj.o
	ar rcs libqg.a $(OBJS)  module_QG_tl.o module_QG_adj.o

QG.exe: QG.f90 $(OBJS)
	$(FC) $(FFLAGS) -o QG.exe QG.f90 $(OBJS) $(INCLUDES) $(LIBS)

QGTLAD.exe: QGTLAD.f90 $(OBJS)
	$(FC) $(FFLAGS) -o QGTLAD.exe QGTLAD.f90 $(OBJS) $(INCLUDES) $(LIBS)

QGTL_test.exe: QGTL_test.f90 $(OBJS)
	$(FC) $(FFLAGS) -o QGTL_test.exe QGTL_test.f90 $(OBJS) $(INCLUDES) $(LIBS)

QGADJ_test.exe: QGADJ_test.f90 $(OBJS)
	$(FC) $(FFLAGS) -o QGADJ_test.exe QGADJ_test.f90 $(OBJS) $(INCLUDES) $(LIBS)

nag.o: nag.f
	$(FC) $(FFLAGS) -c nag.f

nagb.o: nagb.f90
	$(FC) $(FFLAGS) -c nagb.f90

module_QG_model.o: module_QG_model.f90 module_QG_config.o module_QG_ggsp.o nag.o
	$(FC) $(FFLAGS) -c module_QG_model.f90 $(INCLUDES)

module_QG_tl.o: module_QG_tl.f90 module_QG_config.o module_QG_ggsp.o module_QG_util.o nag.o 
	$(FC) $(FFLAGS) -c module_QG_tl.f90 $(INCLUDES)

module_QG_adj.o: module_QG_adj.f90 module_QG_config.o module_QG_ggsp.o module_QG_util.o nag.o 
	$(FC) $(FFLAGS) -convert native -c module_QG_adj.f90 $(INCLUDES)

module_QG_config.o: module_QG_config.f90
	$(FC) $(FFLAGS) -c module_QG_config.f90 $(INCLUDES)

module_QG_writer.o: module_QG_writer.f90 module_QG_model.o
	$(FC) $(FFLAGS) -c module_QG_writer.f90 $(INCLUDES)

module_QGTL_writer.o: module_QGTL_writer.f90 module_QG_tl.o
	$(FC) $(FFLAGS) -c module_QGTL_writer.f90 $(INCLUDES)

module_QGADJ_writer.o: module_QGADJ_writer.f90 module_QG_adj.o
	$(FC) $(FFLAGS) -c module_QGADJ_writer.f90 $(INCLUDES)

module_QG_reader.o: module_QG_reader.f90 module_QG_model.o
	$(FC) $(FFLAGS) -c module_QG_reader.f90 $(INCLUDES)

#module_QG_ggsp.o: module_QG_ggsp.f90 nag.o nagb.o
module_QG_ggsp.o: module_QG_ggsp.f90 nag.o
	$(FC) $(FFLAGS) -c module_QG_ggsp.f90 $(INCLUDES)

module_QG_util.o: module_QG_util.f90
	$(FC) $(FFLAGS) -c module_QG_util.f90 $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod *.a

