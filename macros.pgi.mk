FC=pgfortran
FCOPT=-O3
FCDEBUG=
FFLAGS=-g $(FCOPT) $(FCDEBUG) -Minstrument=functions
LAPACK=-llapack -lblas
GPTLFLAGS=-Minstrument=functions 
GPTL=/contrib/gptl/gptl-v5.5_nompi_noomp_pgi
NETCDF_INCLUDES=-I$(NETCDF)/include
NETCDF_LIBS=-L$(NETCDF)/lib -lnetcdf
