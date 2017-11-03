FC=ifort
FCOPT=-O2 -fp-model precise
FCDEBUG=#-O0 -check all,noarg_temp_created -fpe0
FFLAGS=-g -traceback $(FCOPT) $(FCDEBUG)
#LAPACK=-mkl=parallel
LAPACK=-mkl=sequential
GPTLFLAGS=-finstrument-functions 
#GPTL=/contrib/gptl/gptl-v5.5_nompi_noomp
GPTL=/contrib/gptl/gptl-v5.5_nompi_omp/
NETCDF_INCLUDES=-I$(NETCDF)/include
NETCDF_LIBS=-L$(NETCDF)/lib -lnetcdf -lnetcdff
OPENMP=-qopenmp -qopt-report-stdout -qopt-report=5 -openmp-report  -override_limits
ENDIAN_FLAG=-convert big_endian
