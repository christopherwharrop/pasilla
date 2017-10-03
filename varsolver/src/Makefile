include ../../macros.mk

UTILITIES=../../utilities

INCLUDES=-I. -I$(GPTL)/include            \
         -I../../utilities/include        \
         -I../../models/abstract/include  \
         -I../../models/sine/include      \
         -I../../models/lorenz96/include  \
         -I../../models/QG/include

LIBS=-L$(GPTL)/lib -lgptl                   \
     -L$(NETCDF)/lib -lnetcdf -lnetcdff     \
     -L../../models/sine/lib -lsine         \
     -L../../models/lorenz96/lib -llorenz96 \
     -L../../models/QG/lib -lqg             \
     -L$(UTILITIES)/lib -lpasilla_util      \
     -L../../models/QG/src/tapenade/ADFirstAidKit -lstack


FFLAGS+=$(OPENMP)
FFLAGS+=$(GPTLFLAGS)

all: varsolver_sine.exe varsolver_l96.exe varsolver_qg.exe verify_qg.exe verify_l96.exe

varsolver_sine.exe: varsolver_sine.f90 module_varsolver_sine.o
	$(FC) -o varsolver_sine.exe varsolver_sine.f90 module_varsolver_sine.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver_sine.o: module_varsolver_sine.f90
	$(FC) -c module_varsolver_sine.f90 $(FFLAGS) $(INCLUDES)

varsolver_l96.exe: varsolver_l96.f90 module_varsolver_l96.o
	$(FC) -o varsolver_l96.exe varsolver_l96.f90 module_varsolver_l96.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver_l96.o: module_varsolver_l96.f90
	$(FC) -c module_varsolver_l96.f90 $(FFLAGS) $(INCLUDES)

varsolver_qg.exe: varsolver_qg.f90 module_varsolver_qg.o
	$(FC) -o varsolver_qg.exe varsolver_qg.f90 module_varsolver_qg.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver_qg.o: module_varsolver_qg.f90 ../../models/QG/lib/libqg.a
	$(FC) -c module_varsolver_qg.f90 $(FFLAGS) $(INCLUDES) $(LAPACK) $(NETCDF_INCLUDES)

verify_qg.exe: verify_qg.f90 module_verify_qg.o
	$(FC) -o verify_qg.exe verify_qg.f90 module_verify_qg.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_verify_qg.o: module_verify_qg.f90 ../../models/QG/lib/libqg.a
	$(FC) -c module_verify_qg.f90 $(FFLAGS) $(INCLUDES) $(LAPACK) $(NETCDF_INCLUDES)

verify_l96.exe: verify_l96.f90 module_verify_l96.o
	$(FC) -o verify_l96.exe verify_l96.f90 module_verify_l96.o $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_verify_l96.o: module_verify_l96.f90 ../../models/lorenz96/lib/liblorenz96.a
	$(FC) -c module_verify_l96.f90 $(FFLAGS) $(INCLUDES) $(LAPACK) $(NETCDF_INCLUDES)

clean:
	rm -f *.exe *.o *.mod

