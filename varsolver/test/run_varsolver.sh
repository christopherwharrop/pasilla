#!/bin/sh -l

# Load modules
module load intel/16.1.150
module load netcdf

# Set OMP environment
export OMP_STACKSIZE=1G

# Set locations of data, namelists, and executables
PASILLA_DIR=/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.tapenade
TEST_DIR=${PASILLA_DIR}/models/lorenz96/test
OBS_DIR=${PASILLA_DIR}/models/lorenz96/obs
NATURE_DIR=${PASILLA_DIR}/models/lorenz96/nature  # Needed for initial analysis file
VARSOLVER_DIR=${PASILLA_DIR}/varsolver
LORENZ96_DIR=${PASILLA_DIR}/models/lorenz96

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Brian.Etherton/bin

# Set experiment location
BASE_DIR=${PASILLA_DIR}/models/lorenz96/test

# Set model parameters
lorenz96_forcing=8.10 
lorenz96_time_step=0.000833
start_offset=10

method=$1
sigma=$2
alpha=$3

echo $method, $sigma, $alpha

one_hour=10                      # There are 10 steps in one "hour"
(( fcst_interval=6*$one_hour ))  # Make a forecast every 6 "hours"

f=40060

# Calculate the times 1 hour ahead and 1 hour back
(( fplus1=$f+$one_hour ))
(( fminus1=$f-$one_hour ))

# Calculate the previous forecast analysis time
(( prev_f=$f - $fcst_interval ))

# Make versions of variables with leading zeros
f7=`printf "%07d" $f`
prev_f7=`printf "%07d" $prev_f`

# Set OMP options for this method
export OMP_NUM_THREADS=1
if [ $method -eq 5 ]; then
  export OMP_NUM_THREADS=3
fi

# Create a work directory for the DA and cd into it
workdir=`printf "%07d/$method/%04.2f_%06.4f/varsolverprd" $f $sigma $alpha`
rm -rf $BASE_DIR/$workdir
mkdir -p $BASE_DIR/$workdir
cd $BASE_DIR/$workdir

# Bring over the background files
tidx=1
for t in $fminus1 $f $fplus1; do
  t7=`printf "%07d" $t`
  bkgdir="$TEST_DIR/$prev_f7/$method/lorenz96prd"
  bkgfile="lorenz96out_${t7}.nc"
  ln -s $bkgdir/$bkgfile bkgin_000000${tidx}.nc
  (( tidx=$tidx+1 ))
done

# Bring in the obs file
obsfile="lorenz96obs_${f7}.txt"
ln -s $OBS_DIR/method_${method}/$obsfile lorenz96obs_${method}.txt

# Bring in the namelist and set it up
cp $VARSOLVER_DIR/parm/varsolver_l96.namelist .
sed -i "s/mthd = [[:digit:]]/mthd = ${method}/" varsolver_l96.namelist
sed -i "s/alph = [[:digit:]\.]*/alph = ${alpha}/" varsolver_l96.namelist
sed -i "s/sigma = [[:digit:]\.]*/sigma = ${sigma}/" varsolver_l96.namelist

# Run the assimilation
cp $VARSOLVER_DIR/exe/varsolver_l96.exe .
./varsolver_l96.exe < varsolver_l96.namelist > varsolver_l96.stdout
