FC=ifort
FCOPT=-O3
FCDEBUG=
FFLAGS=-g -traceback $(FCOPT) $(FCDEBUG)
LAPACK=-mkl=sequential
GPTLFLAGS=-finstrument-functions 
GPTL=/contrib/gptl/gptl-v5.5_nompi_noomp
