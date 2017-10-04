#!/bin/sh -l

# Load modules
module load intel/16.1.150
module load netcdf

# Set OMP environment
export OMP_STACKSIZE=1G

# Set locations of data, namelists, and executables
PASILLA_DIR=/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.tapenade
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


one_hour=10                      # There are 10 steps in one "hour"
(( fcst_length=120*$one_hour ))   # Forecast lead time = 120 "hours"
(( fcst_interval=6*$one_hour ))  # Make a forecast every 6 "hours"
start_fcst=40000                 # Starting forecast step
(( end_fcst=$start_fcst+120*24*$one_hour ))  # Perform forecasts for 120 days
#(( end_fcst=$start_fcst ))  # Perform forecasts for 120 days

# Loop over forecasts to perform
f=$start_fcst
while [ $f -le $end_fcst ]; do

  echo "Running forecast for step: $f"

  # Calculate the times 1 hour ahead and 1 hour back
  (( fplus1=$f+$one_hour ))
  (( fminus1=$f-$one_hour ))

  # Calculate the previous forecast analysis time
  (( prev_f=$f - $fcst_interval ))

  # Make versions of variables with leading zeros
  f7=`printf "%07d" $f`
  prev_f7=`printf "%07d" $prev_f`

  # Loop over methods
# for method in 1 2 3 4 5; do
# for method in 1 2 3 4; do
#  for method in 1 2 3; do
  for method in 3; do

    # Set OMP options for this method
    method_dir=$method
    export OMP_NUM_THREADS=1
    if [ $method -eq 5 ]; then
      export OMP_NUM_THREADS=3
    fi

    # If this is not the first forecast, do the DA
    if [ $f -gt $start_fcst ]; then

      # Create a work directory for the DA and cd into it
      workdir=`printf "%07d/$method_dir/varsolverprd" $f`
      mkdir -p $BASE_DIR/$workdir
      cd $BASE_DIR/$workdir

      # Bring over the background files
      tidx=1
      for t in $fminus1 $f $fplus1; do
        t7=`printf "%07d" $t`
        bkgdir="../../../$prev_f7/$method_dir/lorenz96prd"
        bkgfile="lorenz96out_${t7}.nc"
#        cp $bkgdir/$bkgfile .
#        ln -s $bkgfile bkgin_000000${tidx}.nc
        ln -s $bkgdir/$bkgfile bkgin_000000${tidx}.nc
        (( tidx=$tidx+1 ))
      done

      # Bring in the obs file
      obsfile="lorenz96obs_${method}_${f7}.txt"
#      cp $OBS_DIR/method_${method}/$obsfile .
#      ln -s $obsfile lorenz96obs_${method}.txt
      ln -s $OBS_DIR/method_${method}/$obsfile lorenz96obs_${method}.txt

      # Bring in the namelist and set it up
      cp $VARSOLVER_DIR/parm/varsolver_l96.namelist .
      sed -i "s/mthd = [[:digit:]]/mthd = ${method}/" varsolver_l96.namelist

      # Run the assimilation
      cp $VARSOLVER_DIR/exe/varsolver_l96.exe .
      ./varsolver_l96.exe < varsolver_l96.namelist > varsolver_l96.stdout

      # Post-process the timing data
      ${GPTL_PATH}/hex2name.pl ./varsolver_l96.exe ./timing.0 > ./timing.varsolver_l96.txt
      ${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver_l96.txt > timing.varsolver_l96.parse.txt

    fi  # if f >= start_fcst

    # Run the lorenz96 model forecast

    # Create a work directory for the forecast and cd into it
    workdir=`printf "%07d/$method_dir/lorenz96prd" $f`
    mkdir -p $BASE_DIR/$workdir
    cd $BASE_DIR/$workdir

    # If this is the first forecast get analysis from nature
    anlfile=`printf "lorenz96out_%07d.nc" $f`      
    if [ $f -eq $start_fcst  ]; then

      # Get the header for the nature file at fcst_start and fix the forcing value
      natureanl=`printf "lorenz96out_%07d.nc" $f`
      ncdump -h $NATURE_DIR/$natureanl | sed "s/:model_forcing = [^[:blank:]]*/:model_forcing = ${lorenz96_forcing}/" | head -n -1 > anlheader.txt

      # Get the data for the nature file from a different time
      (( foffset=$f+$start_offset ))
      natureanl=`printf "lorenz96out_%07d.nc" $foffset`
      ncdump -v Coordinates,Location,State $NATURE_DIR/$natureanl | grep -A 10000 data: > anldata.txt

      # Put the header and data together into our analysis file
      cat anlheader.txt anldata.txt > anal.txt
      ncgen -o $anlfile anal.txt

    else  # Get analysis from DA
      if [ $method -lt 3 ]; then
        cp ../varsolverprd/bkgout_0000001.nc $anlfile
      else
        cp ../varsolverprd/bkgout_0000002.nc $anlfile
      fi
    fi
  
    # Copy the namelist and set it up
    cp $LORENZ96_DIR/parm/lorenz96.namelist .    
    sed -i "s/forcing = [^[:blank:]]*,/forcing = ${lorenz96_forcing},/" lorenz96.namelist
    sed -i "s/time_step = [^[:blank:]]*,/time_step = ${lorenz96_time_step},/" lorenz96.namelist
    sed -i "s/start_step = [^[:blank:]]*,/start_step = ${f},/" lorenz96.namelist
    sed -i "s/run_steps = [^[:blank:]]*,/run_steps = ${fcst_length},/" lorenz96.namelist
    sed -i "s/output_interval_steps = [^[:blank:]]*,/output_interval_steps = 10,/" lorenz96.namelist

    # Run the lorenz96 model
    cp $LORENZ96_DIR/exe/lorenz96.exe .
    ./lorenz96.exe < ./lorenz96.namelist 2>&1 lorenz96.out

  done # for method

  # Testing/debugging exit
#  if [ $f -gt $start_fcst ]; then
#    exit
#  fi

  # Increment forecast analysis time
  (( f = $f + $fcst_interval ))

done  # while f
