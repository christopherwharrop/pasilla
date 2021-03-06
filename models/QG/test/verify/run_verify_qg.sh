#!/bin/sh -l

PASILLA_DIR=/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.tapenade
NATURE_DIR=${PASILLA_DIR}/models/QG/nature  # Needed for initial analysis file
VARSOLVER_DIR=${PASILLA_DIR}/varsolver
QG_DIR=${PASILLA_DIR}/models/QG

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=${HOME}/bin

# Set experiment location
#BASE_DIR=${QG_DIR}/test/sigma_${parm1}_${parm2}
#BASE_DIR=${QG_DIR}/test/sigma_${sigma}

module load intel/18.1.163
module load netcdf

one_hour=3                        # There are 3 steps in one "hour"
(( fcst_length=120*$one_hour ))   # Forecast lead time = 120 "hours"
(( fcst_interval=12*$one_hour ))  # Make a forecast every 12 "hours"
start_fcst=0                      # Starting forecast step
(( end_fcst=$start_fcst+20*24*$one_hour ))  # Verify forecasts until 20th day

# Link to the input for QG model instantiation
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT106.dat
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT106.dat
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgforcingT106.nc
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qginitT106.nc

ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT21.dat
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT21.dat
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgforcingT21.nc
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qginitT21.nc

h=$1

# Start verifying after 10 days
(( i=$start_fcst + 20*$fcst_interval ))
while [ $i -le $end_fcst ]; do

#  for m in 1 3; do
  for m in 5; do

    v=`expr $i + $h`
    echo $i $h $v

    hh=` printf %07d $h `
    ii=` printf %07d $i `
    vv=` printf %07d $v `
    echo $ii  $hh $vv 

    rm -f truth.nc
#    ln -s /home/Brian.Etherton/da-test/pasilla-QGTEST/models/QG/run/T106/qgout_${vv}.nc truth.nc
    ln -s ${NATURE_DIR}/qgout_${vv}.nc truth.nc

    rm -f model.nc
    ln -s ../../${ii}/${m}/qgprd_assim_off/qgout_${vv}.nc model.nc
    ${VARSOLVER_DIR}/exe/verify_qg.exe > verify_qg.off.stdout 2> verify_qg.off.stderr
    mv verify.txt verify.off.method${m}.init${ii}.lead${hh}.txt

    rm -f model.nc
    ln -s ../../${ii}/${m}/qgprd_assim_on/qgout_${vv}.nc model.nc
    ${VARSOLVER_DIR}/exe/verify_qg.exe  > verify_qg.on.stdout 2> verify_qg.on.stderr
    mv verify.txt verify.on.method${m}.init${ii}.lead${hh}.txt

  done

  (( i= $i + $fcst_interval ))

done 
