include ../../macros.mk

INCLUDES=-I$(GPTL)/include -I../../utilities/include -I../../models/abstract/include -I../../models/sine/include -I../../models/lorenz96/include -I.
LIBS=-L../../models/sine/lib -lsine -L$(NETCDF)/lib -lnetcdf -lnetcdff -L../../models/lorenz96/lib -llorenz96
FFLAGS+=$(OPENMP)
#FFLAGS+=$(GPTLFLAGS)

OBJS=module_config.o \
     module_background.o \
     module_observations.o \
     module_background_covariance.o \
     module_observation_covariance.o \
     module_innovation_vector.o \
     module_observation_operator.o \
     module_solver.o 

varsolver_l96.exe: varsolver.f90 $(OBJS)
	$(FC) -o varsolver_l96.exe varsolver.f90 $(OBJS) $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_config.o: module_config.f90
	$(FC) -c module_config.f90 $(FFLAGS) $(INCLUDES)

module_background.o: module_background.f90
	$(FC) -c module_background.f90 $(FFLAGS) $(INCLUDES)

module_observations.o: module_observations.f90
	$(FC) -c module_observations.f90 $(FFLAGS) $(INCLUDES)

module_background_covariance.o: module_background_covariance.f90
	$(FC) -c module_background_covariance.f90 $(FFLAGS) $(INCLUDES)

module_observation_covariance.o: module_observation_covariance.f90
	$(FC) -c module_observation_covariance.f90 $(FFLAGS) $(INCLUDES)

module_innovation_vector.o: module_innovation_vector.f90 module_background.o module_observations.o module_observation_operator.o
	$(FC) -c module_innovation_vector.f90 $(FFLAGS) $(INCLUDES)

module_observation_operator.o: module_observation_operator.f90
	$(FC) -c module_observation_operator.f90 $(FFLAGS) $(INCLUDES)

module_solver.o: module_solver.f90
	$(FC) -c module_solver.f90 $(FFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod

