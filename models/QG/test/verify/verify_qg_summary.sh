#!/bin/sh -l

h=0
until [ $h -gt 360 ]
do
hh=` printf %03d $h `

m=1
until [ $m -gt 5 ]
do

grep "MSE ALL" ./vx${hh}/v*on*method${m}*init0000*txt | tail -111 > vx${hh}_${m}_on.txt
grep "MSE ALL" ./vx${hh}/v*of*method${m}*init0000*txt | tail -111 > vx${hh}_${m}_off.txt

m=`expr $m + 1 `
done

h=`expr $h + 36 ` 
done
