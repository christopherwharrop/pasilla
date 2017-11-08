#!/bin/sh -l

module load intel/16.1.150
module load netcdf
module load grads

# Create psidiff.z[123].dat files for input to bmatrix.exe
if [ ! -f psidiff.z1.dat ]; then
  grads -blc plot.grads.z1.gs &
fi
if [ ! -f psidiff.z2.dat ]; then
  grads -blc plot.grads.z2.gs &
fi
if [ ! -f psidiff.z3.dat ]; then
  grads -blc plot.grads.z3.gs &
fi

# Wait for background Grads scripts to complete
wait

# Make the latlon.dat file if needed
if [ ! -f latlon.dat ]; then
  ifort -o latlon.exe latlon.f90 -O3
  ./latlon.exe
fi

# Make B squared
ifort -o bmatrix.exe bmatrix.f90 -O3 -mkl=parallel -I${NETCDF}/include -L${NETCDF}/lib -lnetcdff -lnetcdf
./bmatrix.exe > bmatrix.out

# Make B
ifort -o sqrt_b.exe sqrt_b.f90 -O3 -mkl=parallel -I${NETCDF}/include -L${NETCDF}/lib -lnetcdff -lnetcdf
./sqrt_b.exe > sqrt_b.out


