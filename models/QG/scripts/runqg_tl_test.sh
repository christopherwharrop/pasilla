#!/bin/bash -l

# This scripts runs the QG model.
#
# resol sets the resolution
# qgdir must point to the directory of the distribution
# expid identifies the directory where the output of the run will be stored
# obsfile is the datafile with observations used to calculate the forcing from

#module load intel/16.1.150
#module load netcdf
#module load pgi/16.10
#module load netcdf/4.4.0

# Set the path to GPTL
GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin

# Set the path to the input data
#INPUT_DATA=/lfs3/projects/gsd-hpcs/QG/inputdata
#INPUT_DATA=$WORK/QG/inputdata
INPUT_DATA=/scratch4/BMC/gsd-hpcs/QG/inputdata

# Set the resolution of the run
resol="21"

# Set the experiment id
expid='tltest'

# Set the path to the QG model install - DO NOT CHANGE THIS LINE
qgdir=$(dirname $(dirname $(readlink -m $0) ))

# Set the path to the runtime parameter files
parmdir="${qgdir}/parm"

# Set the path to the run directory
rundir="${qgdir}/run/${expid}"

# Set the name of the obs file to use - This must always be the T106 file (for now)
obsfile="sf7910T106.shfs"

# Create an empty run directory and cd into it
rm -rf ${rundir}
mkdir -p ${rundir}
cd ${rundir}

# Copy the inputdata into place
cp -prd ${INPUT_DATA}/${obsfile} .
cp -prd ${INPUT_DATA}/qgcoefT${resol}.dat .
cp -prd ${INPUT_DATA}/qgbergT${resol}.dat .

# These are needed for boostrapping model state without a restart file
cp -prd ${INPUT_DATA}/qginitT${resol}.nc .
cp -prd ${INPUT_DATA}/qgforcingT${resol}.nc .

# Make sure obs file exists.
if [ ! -e $obsfile ]; then
  echo "${qgdir}/inputdata/$obsfile does not exist" ; exit 1
fi

# Copy the namelist into the run directory
cp ${parmdir}/QG.namelist QG.namelist

# Set the resolution in the namelist
sed -i "s/resolution = [[:digit:]]*/resolution = ${resol}/" QG.namelist

# Set the obs file in the namelist
sed -i "s/obsfile = '.*'/obsfile = \'${obsfile}\'/" QG.namelist

# Set the run steps in the namelist
sed -i "s/run_steps = [[:digit:]]*/run_steps = 8/" QG.namelist

# Copy the executable to the run directory
cp ${qgdir}/exe/QGTL_test.exe .

# Run the model
./QGTL_test.exe < QG.namelist

# Parse the GPTL timing
#${GPTL_PATH}/hex2name.pl ./QGTLAD.exe ./timing.0 > ./timing.QG.txt
