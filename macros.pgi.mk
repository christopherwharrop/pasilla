FC=pgfortran
FCOPT=-O2 -Kieee
FCDEBUG=#-O0 -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -traceback
FFLAGS=-g $(FCOPT) $(FCDEBUG)
LAPACK=-llapack -lblas
GPTLFLAGS=-Minstrument=functions 
GPTL=/contrib/gptl/gptl-v5.5_nompi_noomp_pgi
NETCDF_INCLUDES=-I$(NETCDF)/include
NETCDF_LIBS=-L$(NETCDF)/lib -lnetcdf -lnetcdff
OPENMP=-mp
