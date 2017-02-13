#!/bin/bash -l

# This scripts runs the QG model.
#
# resol sets the resolution
# qgdir must point to the directory of the distribution
# expid identifies the directory where the output of the run will be stored
# obsfile is the datafile with observations used to calculate the forcing from

module load intel/16.1.150
module load netcdf
#module load pgi/16.10
#module load netcdf/4.4.0

# Set the resolution of the run
resol="21"

# Set the experiment id
expid='harrop'

# Set base paths
qgdir="/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.top/models/QG"
parmdir="${qgdir}/parm"
rundir="${qgdir}/run/${expid}"
obsfile="sf7910T106.shfs"

# Create an empty run directory and cd into it
rm -rf ${rundir}
mkdir -p ${rundir}
cd ${rundir}

# Copy the inputdata into place
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/${obsfile} .
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT${resol}.dat .
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT${resol}.dat .
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/qgstartT${resol}.dat .

# Make sure obs file exists.
if [ ! -e $obsfile ]; then
  echo "${qgdir}/inputdata/$obsfile does not exist" ; exit 1
fi

# Copy the namelist into the run directory
cp ${parmdir}/namelist namelist.input

# Set the resolution in the namelist
sed -i "s/resolution = [[:digit:]]*/resolution = ${resol}/" namelist.input

# Set the obs file in the namelist
sed -i "s/obsfile = '.*'/obsfile = \'${obsfile}\'/" namelist.input

# Copy the executable to the run directory
cp ${qgdir}/exe/QG.exe .

# Run the model
./QG.exe
