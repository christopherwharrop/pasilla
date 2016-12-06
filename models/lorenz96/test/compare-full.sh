#!/bin/sh

t=0
while [ $t -lt 24005 ]
do

for c in '1' '2' '3' '4_1' '4_3'  
do

i=0
while [ $i -lt 785 ]
#while [ $i -lt 15 ]
do
vt=$(printf %06d `expr 40000 + $t + $i`)
it=$(printf %06d `expr 40000 + $t`)
echo -ne ${vt} '\t' ${it} '\t' ${i} '\t' ${c} '\t' 
#echo "./compare2.sh nature/lorenz96out_0${vt}.nc 0${it}/${c}/lorenz96prd/lorenz96out_0${vt}.nc comp.nc"
./compare2.sh nature/lorenz96out_0${vt}.nc 0${it}/${c}/lorenz96prd/lorenz96out_0${vt}.nc comp.nc 
i=`expr $i + 60`
done

done

t=`expr $t + 60`
done

