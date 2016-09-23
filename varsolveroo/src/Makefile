include ../../macros.mk

INCLUDES=-I$(GPTL)/include -I../../utilities/include -I.
LIBS=-L$(GPTL)/lib -lgptl
FFLAGS+=-openmp

OBJS=module_varsolver.o module_background.o module_observations.o module_background_covariance.o module_observation_covariance.o

varsolver.exe: varsolver.f90 $(OBJS)
	$(FC) -o varsolver.exe varsolver.f90 $(OBJS) $(FFLAGS) $(LAPACK) $(INCLUDES) $(LIBS)

module_varsolver.o: module_varsolver.f90 module_background.o module_observations.o module_background_covariance.o module_observation_covariance.o
	$(FC) -c module_varsolver.f90 $(FFLAGS) $(GPTLFLAGS) $(INCLUDES)

module_background.o: module_background.f90
	$(FC) -c module_background.f90 $(FFLAGS) $(GPTLFLAGS) $(INCLUDES)

module_observations.o: module_observations.f90
	$(FC) -c module_observations.f90 $(FFLAGS) $(GPTLFLAGS) $(INCLUDES)

module_background_covariance.o: module_background_covariance.f90
	$(FC) -c module_background_covariance.f90 $(FFLAGS) $(GPTLFLAGS) $(INCLUDES)

module_observation_covariance.o: module_observation_covariance.f90
	$(FC) -c module_observation_covariance.f90 $(FFLAGS) $(GPTLFLAGS) $(INCLUDES)

clean:
	rm -f *.exe *.o *.mod

