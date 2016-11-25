#!/bin/sh
rm comp.nc
rm state.nc

#module load intel/16.1.150 netcdf nco

ncdiff -v State $1 $2 comp.nc
ncwa -O -y rms comp.nc state.nc
ncdump state.nc | grep "State =" | cut -d" " -f4
