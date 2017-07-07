#!/bin/sh
rm -f comp.nc
rm -f state.nc

#module load intel/16.1.150 netcdf nco

ncdiff -v Psig $1 $2 comp.nc
ncwa -O -y rms comp.nc state.nc
ncdump state.nc | grep "Psig =" | cut -d" " -f4
