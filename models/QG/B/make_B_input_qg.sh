#!/bin/sh -l

module load intel/16.1.150
module load netcdf
module load grads

# Remove existing input files
rm -f psidiff.z*.dat latlon.dat

# Create psidiff.z[123].dat files for input to bmatrix.exe
echo "Creating psidiff.z1.dat"
grads -blc plot.grads.z1.gs &
echo "Creating psidiff.z2.dat"
grads -blc plot.grads.z2.gs &
echo "Creating psidiff.z3.dat"
grads -blc plot.grads.z3.gs &

# Wait for background Grads scripts to complete
wait

# Make the latlon.dat file if needed
echo "Creating latlon.dat"
./latlon.exe
