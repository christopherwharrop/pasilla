#!/bin/sh -l

# Get sigma for the experiment to verify
sigma=$1
#parm1=$1
#parm2=$2

PASILLA_DIR=/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.tapenade
NATURE_DIR=${PASILLA_DIR}/models/QG/nature  # Needed for initial analysis file
VARSOLVER_DIR=${PASILLA_DIR}/varsolver
QG_DIR=${PASILLA_DIR}/models/QG

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=${HOME}/bin

# Set experiment verification location
BASE_DIR=${QG_DIR}/test/sigma_$sigma/verify
#BASE_DIR=${QG_DIR}/test/sigma_${parm1}_${parm2}/verify

# Create experiment verification location
rm -rf $BASE_DIR
mkdir $BASE_DIR
cd $BASE_DIR

# Set up run parameters
one_hour=3                       # There are 3 steps in one "hour"
(( fcst_length=120*$one_hour ))  # Forecast lead time = 120 "hours"
(( fcst_interval=12*$one_hour )) # Make a forecast every 12 "hours"

# Verify each lead time
h=0
while [ $h -le $fcst_length ]; do

  hh=` printf %04d $h `

  rm -rf vx${hh}
  mkdir -p vx${hh}
  cd vx${hh}
  cp ${QG_DIR}/test/verify/run_verify_qg.sh .
  ./run_verify_qg.sh ${h} > verify_qg_lead${hh}.log &
  cd ../

  (( h = $h + $fcst_interval ))

done

# Wait for verification of each lead time to complete
wait

# Run summary scripts
cd $BASE_DIR
${QG_DIR}/test/verify/verify_qg_summary.sh
${QG_DIR}/test/verify/verify_qg_summary.py
