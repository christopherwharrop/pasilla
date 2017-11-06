FC=ifort
FCOPT=-O2 -fp-model precise
FCDEBUG=#-O0 -check all,noarg_temp_created -fpe0 -init=snan,arrays
LAPACK=-mkl=parallel
FFLAGS=-g -traceback $(FCOPT) $(LAPACK) $(FCDEBUG)
GPTLFLAGS=-finstrument-functions 
GPTL=/contrib/gptl/gptl-v5.5_nompi_omp/
NETCDF_INCLUDES=-I$(NETCDF)/include
NETCDF_LIBS=-L$(NETCDF)/lib -lnetcdf -lnetcdff
OPENMP=-qopenmp -opt-report-file=stdout -opt-report-phase:openmp
ENDIAN_FLAG=-convert big_endian
