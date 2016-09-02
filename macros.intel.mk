FC=ifort
FCOPT=-O3
FCDEBUG=
FFLAGS=-g -traceback $(FCOPT) $(FCDEBUG)
LAPACK=-mkl=parallel
GPTLFLAGS=-finstrument-functions 
GPTL=/contrib/gptl/gptl-v5.5_nompi_noomp
NETCDF_INCLUDES=-I$(NETCDF)/include
NETCDF_LIBS=-L$(NETCDF)/lib -lnetcdf -lnetcdff
