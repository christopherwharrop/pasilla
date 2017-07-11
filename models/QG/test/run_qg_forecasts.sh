#!/bin/sh -l

# Load modules
module load intel/16.1.150
module load netcdf

# Check usage
if [ $# -ne 2 ]; then
  echo "USAGE: run_qg_forecasts.sh START_FCST END_FCST"
  exit 1
else
  start_fcst=$1                    # Starting forecast step
  end_fcst=$2                      # Perform forecasts for 120 days
fi

# Set OMP environment
export OMP_STACKSIZE=1G

# Set locations of data, namelists, and executables
PASILLA_DIR=$(dirname $(dirname $(dirname $(dirname $(readlink -m $0)))))  # DO NOT CHANGE THIS
OBS_DIR=${PASILLA_DIR}/models/QG/obs
NATURE_DIR=${PASILLA_DIR}/models/QG/nature  # Needed for initial analysis file
VARSOLVER_DIR=${PASILLA_DIR}/varsolver
QG_DIR=${PASILLA_DIR}/models/QG

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=${HOME}/bin

# Set experiment location
BASE_DIR=${QG_DIR}/test

# Set model parameters
start_step=0
spinup_steps=0
run_steps=360
output_interval_steps=18
resolution=21
time_step=1200

one_hour=3                       # There are 3 steps in one "hour"
(( fcst_length=120*$one_hour ))  # Forecast lead time = 78 "hours"
(( fcst_interval=6*$one_hour ))  # Make a forecast every 6 "hours"
start_epoch=0                    # First forecast in the series (cold start cycle)
(( assimilation_window=12*$one_hour ))       # Assimilate observations every 12 "hours"

# Loop over forecasts to perform
f=$start_fcst
while [ $f -le $end_fcst ]; do

  echo "Running forecast for step: $f"

  # Calculate the times 1 hour ahead and 1 hour back
  (( fplus1=$f+$assimilation_window ))
  (( fminus1=$f-$assimilation_window ))

  # Calculate the previous forecast analysis time
  (( prev_f=$f - $assimilation_window ))

  # Make versions of variables with leading zeros
  f7=`printf "%07d" $f`
  prev_f7=`printf "%07d" $prev_f`

  # Loop over methods
# for method in 1 2 3 4 5; do
  for method in 1 3 4 5; do

    # Set OMP options for this method
    method_dir=$method
    export OMP_NUM_THREADS=4
    if [ $method -eq 5 ]; then
      export OMP_NUM_THREADS=12
    fi

    # If this is not the first forecast, do the DA
    if [ $f -gt $start_epoch ]; then

      # Create a work directory for the DA and cd into it
      workdir=`printf "%07d/$method_dir/varsolverprd" $f`
      mkdir -p $BASE_DIR/$workdir
      cd $BASE_DIR/$workdir

      # Bring over the background files
      tidx=1
      for t in $fminus1 $f $fplus1; do
        t7=`printf "%07d" $t`
        bkgdir="../../../$prev_f7/$method_dir/qgprd_assim_on"
        bkgfile="qgout_${t7}.nc"
        ln -s $bkgdir/$bkgfile bkgin_000000${tidx}.nc
        (( tidx=$tidx+1 ))
      done

      # Bring in the obs file BUG!! NEED DIFFERENT OBS FILES FOR EACH METHOD
      obsfile="qgobs_${f7}.txt"
      ln -s $OBS_DIR/$obsfile qgobs_${method}.txt

      # Bring in the namelist and set it up
      cp $VARSOLVER_DIR/parm/varsolver_qg.namelist .
      sed -i "s/mthd = [[:digit:]]/mthd = ${method}/" varsolver_qg.namelist

      # Copy the inputdata into place
      obsfile="sf7910T106.shfs"
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qginitT${resolution}.nc
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgforcingT${resolution}.nc
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/${obsfile}
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT${resolution}.dat
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT${resolution}.dat
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/bT${resolution}.nc b.nc
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/mu_operator.dat

      # Run the assimilation
      cp $VARSOLVER_DIR/exe/varsolver_qg.exe .
      ./varsolver_qg.exe < varsolver_qg.namelist > varsolver_qg.stdout

      # Post-process the timing data
      ${GPTL_PATH}/hex2name.pl ./varsolver_qg.exe ./timing.0 > ./timing.varsolver_qg.txt
      ${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver_qg.txt > timing.varsolver_qg.parse.txt

    fi  # if f >= start_epoch

    # Run the QG model forecast with assimilation

    # Create a work directory for the forecast and cd into it
    workdir=`printf "%07d/$method_dir/qgprd_assim_on" $f`
    mkdir -p $BASE_DIR/$workdir
    cd $BASE_DIR/$workdir

    # Set name of analysis file
    anlfile=`printf "qgout_%07d.nc" $f`      

    # If this is the first forecast get analysis from nature
    if [ $f -eq $start_epoch ]; then

      # These are needed for boostrapping model state without a restart file
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qginitT${resolution}.nc
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgforcingT${resolution}.nc

    else  # Get analysis from DA
      if [ $method -lt 3 ]; then
        cp ../varsolverprd/bkgout_0000001.nc $anlfile
      else
        cp ../varsolverprd/bkgout_0000002.nc $anlfile
      fi
    fi
  
    # Copy the namelist and set it up
    cp $QG_DIR/parm/QG.namelist .    
    sed -i "s/start_step = [^[:blank:]]*,/start_step = ${f},/" QG.namelist
    sed -i "s/spinup_steps = [^[:blank:]]*,/spinup_steps = ${spinup_steps},/" QG.namelist
    sed -i "s/run_steps = [^[:blank:]]*,/run_steps = ${fcst_length},/" QG.namelist
    sed -i "s/output_interval_steps = [^[:blank:]]*,/output_interval_steps = ${output_interval_steps},/" QG.namelist
    sed -i "s/resolution = [^[:blank:]]*,/resolution = ${resolution}/" QG.namelist
    sed -i "s/time_step = [^[:blank:]]*,/time_step = ${time_step}/" QG.namelist

    # Copy the inputdata into place
    obsfile="sf7910T106.shfs"
    ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/${obsfile}
    ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT${resolution}.dat
    ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT${resolution}.dat

    # Run the QG model
    cp $QG_DIR/exe/QG.exe .
    ./QG.exe < ./QG.namelist 

    # Run the QG model forecast withOUT assimilation

    # Create a work directory for the forecast and cd into it
    workdir=`printf "%07d/$method_dir/qgprd_assim_off" $f`
    mkdir -p $BASE_DIR/$workdir
    cd $BASE_DIR/$workdir

    # If this is the first forecast get analysis from nature
    if [ $f -eq $start_epoch  ]; then

      # These are needed for boostrapping model state without a restart file
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qginitT${resolution}.nc
      ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgforcingT${resolution}.nc

    else  # Get analysis from non-assimilated forecast
      cp ../varsolverprd/bkgin_0000002.nc $anlfile
    fi
  
    # Copy the namelist and set it up
    cp $QG_DIR/parm/QG.namelist .    
    sed -i "s/start_step = [^[:blank:]]*,/start_step = ${f},/" QG.namelist
    sed -i "s/spinup_steps = [^[:blank:]]*,/spinup_steps = ${spinup_steps},/" QG.namelist
    sed -i "s/run_steps = [^[:blank:]]*,/run_steps = ${fcst_length},/" QG.namelist
    sed -i "s/output_interval_steps = [^[:blank:]]*,/output_interval_steps = ${output_interval_steps},/" QG.namelist
    sed -i "s/resolution = [^[:blank:]]*,/resolution = ${resolution}/" QG.namelist
    sed -i "s/time_step = [^[:blank:]]*,/time_step = ${time_step}/" QG.namelist

    # Copy the inputdata into place
    obsfile="sf7910T106.shfs"
    ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/${obsfile}
    ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT${resolution}.dat
    ln -s /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT${resolution}.dat

    # Run the QG model
    cp $QG_DIR/exe/QG.exe .
    ./QG.exe < ./QG.namelist 

  done # for method

  # Testing/debugging exit
#  if [ $f -gt $start_fcst ]; then
#    exit
#  fi

  # Increment forecast analysis time
#  (( f = $f + $fcst_interval ))
  (( f = $f + $assimilation_window ))

done  # while f
