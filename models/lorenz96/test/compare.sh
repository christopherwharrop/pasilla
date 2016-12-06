#!/bin/sh

t=1020
while [ $t -lt 4800 ]
do
for c in '1' '2' '3' '4_1' '4_3'  
do
tt=`expr $t + 0`
echo -ne ${t} '\t' ${tt} '\t' ${c} '\t' 
./compare2.sh nature/lorenz96out_004${t}.nc 004${t}/${c}/lorenz96prd/lorenz96out_004${t}.nc comp.nc
done
t=`expr $t + 120`
done

t=1020
while [ $t -lt 4800 ]
do
for c in '1' '2' '3' '4_1' '4_3'
do
tt=`expr $t + 60`
echo -ne ${t} '\t' ${tt} '\t' ${c} '\t' 
./compare2.sh nature/lorenz96out_004${tt}.nc 004${t}/${c}/lorenz96prd/lorenz96out_004${tt}.nc comp.nc
done
t=`expr $t + 120`
done
