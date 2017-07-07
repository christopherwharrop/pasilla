#!/bin/sh

#t=1020
t=0
while [ $t -le 72 ]; do
  for c in '1' '2' '3' '4_1' '4_3'; do
    tt=`expr $t + 0`
    t7=`printf "%07d" $t`
    tt7=`printf "%07d" $tt`
    echo -ne ${t7} '\t' ${tt7} '\t' ${c} '\t' 
    #./compare2.sh nature/lorenz96out_004${t}.nc 004${t}/${c}/lorenz96prd/lorenz96out_004${t}.nc comp.nc
    ./compare2.sh ../nature/qgout_${t7}.nc ${t7}/${c}/qgprd_assim_on/qgout_${t7}.nc
  done
  t=`expr $t + 36`
done

t=1020
while [ $t -le 72 ]; do
  for c in '1' '2' '3' '4_1' '4_3'; do
    tt=`expr $t + 60`
    echo -ne ${t} '\t' ${tt} '\t' ${c} '\t' 
    #./compare2.sh nature/lorenz96out_004${tt}.nc 004${t}/${c}/lorenz96prd/lorenz96out_004${tt}.nc comp.nc
    ./compare2.sh ../nature/qgout_${tt}.nc ${t}/${c}/qgprd/qgout_${tt}.nc
  done
  t=`expr $t + 36`
done
