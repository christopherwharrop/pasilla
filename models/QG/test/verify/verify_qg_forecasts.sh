#!/bin/sh -l

h=0
until [ $h -gt 360 ]
do
hh=` printf %03d $h `

mkdir ./vx${hh}
cd ./vx${hh}
cp ../run_verify_qg.sh .
./run_verify_qg.sh ${h} > verify_qg_lead${hh}.log &
cd ../

h=`expr $h + 36 ` 
done
