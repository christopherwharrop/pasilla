#!/bin/sh -l

module load intel/16.1.150
module load netcdf
module load grads

# Get base directory
BASE_DIR=/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.tapenade/models/QG/B

# Get desired B matrix sigma from command line
#sigma=$1
parm1=$1
parm2=$2

# Create work directory and cd into it
#workdir=${BASE_DIR}/sigma_${sigma}
workdir=${BASE_DIR}/sigma_${parm1}_${parm2}
rm -rf ${workdir}
mkdir -p ${workdir}
cd ${workdir}

# link to input files
ln -s ../psidiff.z1.dat
ln -s ../psidiff.z2.dat
ln -s ../psidiff.z3.dat
ln -s ../latlon.dat

# Copy namelist and modify it
cp ../bmatrix.namelist .
#sed -i "s/sigma = [^[:blank:]]*,/sigma = ${sigma},/" bmatrix.namelist
sed -i "s/parm1 = [^[:blank:]]*,/parm1 = ${parm1},/" bmatrix.namelist
sed -i "s/parm2 = [^[:blank:]]*,/parm2 = ${parm2},/" bmatrix.namelist

# Make B squared
rm -f b_sqr.nc
ln -s ../bmatrix.exe
echo "Running bmatrix.exe for sigma=$sigma"
./bmatrix.exe < bmatrix.namelist > bmatrix.out

# Make B
rm -f b.nc
ln -s ../sqrt_b.exe
echo "Running sqrt_b.exe"
./sqrt_b.exe > sqrt_b.out

# Copy B to the input data directory
#cp b.nc /scratch4/BMC/gsd-hpcs/QG/inputdata/bT21.${sigma}.nc
cp b.nc /scratch4/BMC/gsd-hpcs/QG/inputdata/bT21.${parm1}_${parm2}.nc
