include ../../../macros.mk

UTILITIES=../../../utilities

INCLUDES=-I$(GPTL)/include -I$(UTILITIES)/include $(NETCDF_INCLUDES) -I../../abstract/include -I.
LIBS=-L$(GPTL)/lib -lgptl -L$(UTILITIES)/lib -lpasilla_util $(NETCDF_LIBS)

OBJS=module_L96_model.o module_L96_tl.o module_L96_adj.o module_L96_config.o module_L96_reader.o module_L96_writer.o module_L96TL_writer.o module_L96ADJ_writer.o module_netcdf_writer.o ../../abstract/src/module_model.o

FFLAGS+=$(GPTLFLAGS)

#all: lorenz96.exe lorenz96TL_test.exe lorenz96ADJ_test.exe harvest_lorenz96_obs.exe liblorenz96.a
all: lorenz96.exe L96TLAD.exe liblorenz96.a lorenz96TL_test.exe lorenz96ADJ_test.exe

liblorenz96.a: $(OBJS) module_L96_tl.o module_L96_adj.o
	ar rcs liblorenz96.a $(OBJS) module_L96_tl.o module_L96_adj.o

harvest_lorenz96_obs.exe: module_lorenz96.o harvest_lorenz96_obs.f90
	$(FC) $(FFLAGS) -o harvest_lorenz96_obs.exe harvest_lorenz96_obs.f90 module_lorenz96.o $(INCLUDES) $(LIBS)

#lorenz96.exe: module_lorenz96.o lorenz96.f90
#	$(FC) $(FFLAGS) -o lorenz96.exe lorenz96.f90 module_lorenz96.o $(INCLUDES) $(LIBS)

lorenz96.exe: L96.f90 $(OBJS)
	$(FC) $(FFLAGS) -o lorenz96.exe L96.f90 $(OBJS) $(INCLUDES) $(LIBS)
module_L96_model.o: module_L96_config.o module_L96_model.f90
	$(FC) $(FFLAGS) -c module_L96_model.f90 $(INCLUDES) $(LIBS)
module_L96_tl.o: module_L96_config.o module_L96_tl.f90
	$(FC) $(FFLAGS) -c module_L96_tl.f90 $(INCLUDES) $(LIBS)
module_L96_adj.o: module_L96_config.o module_L96_adj.f90
	$(FC) $(FFLAGS) -c module_L96_adj.f90 $(INCLUDES) $(LIBS)
module_L96_config.o: module_L96_config.f90
	$(FC) $(FFLAGS) -c module_L96_config.f90 $(INCLUDES) $(LIBS)
module_L96_reader.o: module_L96_reader.f90
	$(FC) $(FFLAGS) -c module_L96_reader.f90 $(INCLUDES) $(LIBS)
module_L96_writer.o: module_L96_writer.f90 module_netcdf_writer.o
	$(FC) $(FFLAGS) -c module_L96_writer.f90 $(INCLUDES) $(LIBS)
module_L96TL_writer.o: module_L96TL_writer.f90 module_netcdf_writer.o module_L96_tl.o
	$(FC) $(FFLAGS) -c module_L96TL_writer.f90 $(INCLUDES) $(LIBS)
module_L96ADJ_writer.o: module_L96ADJ_writer.f90 module_netcdf_writer.o module_L96_adj.o
	$(FC) $(FFLAGS) -c module_L96ADJ_writer.f90 $(INCLUDES) $(LIBS)
module_netcdf_writer.o: module_netcdf_writer.f90
	$(FC) $(FFLAGS) -c module_netcdf_writer.f90 $(INCLUDES) $(LIBS)

L96TLAD.exe: L96TLAD.f90 $(OBJS)
	$(FC) $(FFLAGS) -o L96TLAD.exe L96TLAD.f90 $(OBJS) $(INCLUDES) $(LIBS)


lorenz96TL_test.exe: lorenz96TL_test.f90 $(OBJS)
	$(FC) $(FFLAGS) -o lorenz96TL_test.exe lorenz96TL_test.f90 $(OBJS) $(INCLUDES) $(LIBS)

lorenz96ADJ_test.exe: lorenz96ADJ_test.f90
	$(FC) $(FFLAGS) -o lorenz96ADJ_test.exe lorenz96ADJ_test.f90 $(OBJS) $(INCLUDES) $(LIBS)

#module_lorenz96.o: module_lorenz96.f90
#	$(FC) $(FFLAGS) -c module_lorenz96.f90 $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod *.a

