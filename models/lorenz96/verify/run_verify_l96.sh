#!/bin/sh -l

module load intel/16.1.150
module load netcdf

one_hour=10                       # There are 10 steps in one "hour"
(( fcst_length=120*$one_hour ))   # Forecast lead time = 120 "hours"
(( fcst_interval=6*$one_hour ))  # Make a forecast every 6 "hours"
start_fcst=40000                  # Starting forecast step
(( end_fcst=$start_fcst+100*24*$one_hour ))  # Perform forecasts for 120 days

h=$1

i=$start_fcst
while [ $i -le $end_fcst ]; do

  for m in 1 3; do

    v=`expr $i + $h`
    echo $i $h $v

    hh=` printf %07d $h `
    ii=` printf %07d $i `
    vv=` printf %07d $v `
    echo $ii  $hh $vv 
    rm -f truth.nc
    #ln -s /home/Brian.Etherton/da-test/pasilla-QGTEST/models/QG/run/T106/qgout_${vv}.nc truth.nc
    ln -s ../../nature/lorenz96out_${vv}.nc truth.nc
    rm -f model.nc
    #ln -s ../../${ii}/${m}/qgprd_assim_off/qgout_${vv}.nc model.nc
    ln -s ../../test/${ii}/${m}/lorenz96prd/lorenz96out_${vv}.nc model.nc
    rm -f verify.txt
    ../../../../varsolver/exe/verify_l96.exe
    ls
    mv verify.txt verify.method${m}.init${ii}.lead${hh}.txt

  done

  (( i= $i + $fcst_interval ))

done 
