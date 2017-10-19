#!/bin/sh -l

one_hour=3                       # There are 3 steps in one "hour"
(( fcst_length=120*$one_hour ))  # Forecast lead time = 120 "hours"
(( fcst_interval=12*$one_hour )) # Make a forecast every 12 "hours"
start_fcst=0                     # Starting forecast step
#(( end_fcst=$start_fcst+100*24*$one_hour ))  # Perform forecasts for 120 days
(( end_fcst=$start_fcst+20*24*$one_hour ))  # Perform forecasts for 120 days


h=0
while [ $h -le $fcst_length ]; do

  hh=` printf %04d $h `

  rm -rf vx${hh}
  mkdir -p vx${hh}
  cd vx${hh}
  cp ../run_verify_qg.sh .
  ./run_verify_qg.sh ${h} > verify_qg_lead${hh}.log &
#  ./run_verify_qg.sh ${h} > verify_qg_lead${hh}.log 
  cd ../

  (( h = $h + $fcst_interval ))

done
