#!/bin/sh -l

one_hour=3                       # There are 3 steps in one "hour"
(( fcst_length=120*$one_hour ))  # Forecast lead time = 120 "hours"
(( fcst_interval=12*$one_hour )) # Make a forecast every 12 "hours"

h=0
while [ $h -le $fcst_length ]; do

  hh=`printf %04d $h`

  for m in 5; do
#  for m in 1 3; do
#  for m in 1 2 3; do
#  for m in 1; do

#    grep "^ RMSE ALL" ./vx${hh}/v*on*method${m}*init0*txt > vx${hh}_${m}_on.txt
#    grep "^ RMSE ALL" ./vx${hh}/v*of*method${m}*init0*txt > vx${hh}_${m}_off.txt
    grep "^ MSE ALL" ./vx${hh}/v*on*method${m}*init0*txt > vx${hh}_${m}_on.txt
    grep "^ MSE ALL" ./vx${hh}/v*of*method${m}*init0*txt > vx${hh}_${m}_off.txt

  done

  (( h = $h + $fcst_interval ))

done
