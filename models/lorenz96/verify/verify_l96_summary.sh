#!/bin/sh -l

one_hour=10                       # There are 10 steps in one "hour"
(( fcst_length=120*$one_hour ))   # Forecast lead time = 120 "hours"
(( fcst_interval=6*$one_hour ))  # Make a forecast every 6 "hours"
start_fcst=40000                  # Starting forecast step
(( end_fcst=$start_fcst+100*24*$one_hour ))  # Perform forecasts for 120 days

h=0
while [ $h -le $fcst_length ]; do

  hh=`printf %04d $h`

  for m in 1 3; do

    grep "MSE ALL" ./vx${hh}/v*method${m}*init0*txt > vx${hh}_${m}.txt
    #grep "MSE ALL" ./vx${hh}/v*on*method${m}*init0000*txt | tail -111 > vx${hh}_${m}_on.txt
    #grep "MSE ALL" ./vx${hh}/v*of*method${m}*init0000*txt | tail -111 > vx${hh}_${m}_off.txt

  done

  (( h = $h + $fcst_interval ))

done
