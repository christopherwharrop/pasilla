FC=ifort
FCOPT=-O2 -fp-model precise # -xCORE-AVX2
FCDEBUG=#-O0 -check all,noarg_temp_created -fpe0 -init=snan,arrays
LAPACK=-mkl=parallel
FFLAGS=-g -traceback $(FCOPT) $(LAPACK) $(FCDEBUG)
GPTLFLAGS=-finstrument-functions 
GPTL=/contrib/gptl/gptl-v5.5_nompi_omp/
NETCDF_INCLUDES=-I$(NETCDF)/include
NETCDF_LIBS=-L$(NETCDF)/lib -lnetcdf -lnetcdff
OPENMP=-qopenmp -qopt-report-file=stdout -qopt-report-phase:openmp
ENDIAN_FLAG=-convert big_endian
