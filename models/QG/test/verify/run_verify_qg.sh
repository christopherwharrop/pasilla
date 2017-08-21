#!/bin/sh -l

module load intel/16.1.150
module load netcdf

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

i=0
until [ $i -gt 2630 ]
do

m=1
until [ $m -gt 5 ]
do
v=`expr $i + $h`
echo $i $h $v

hh=` printf %07d $h `
ii=` printf %07d $i `
vv=` printf %07d $v `
echo $ii  $hh $vv 
rm truth.nc
ln -s /home/Brian.Etherton/da-test/pasilla-QGTEST/models/QG/run/T106/qgout_${vv}.nc truth.nc
rm model.nc
ln -s ../../${ii}/${m}/qgprd_assim_off/qgout_${vv}.nc model.nc
./verify_qg.exe  > verify_qg.on.stdout 2> verify_qg.on.stderr
mv verify.txt verify.off.method${m}.init${ii}.lead${hh}.txt

rm model.nc
ln -s ../../${ii}/${m}/qgprd_assim_on/qgout_${vv}.nc model.nc
./verify_qg.exe  > verify_qg.off.stdout 2> verify_qg.off.stderr
mv verify.txt verify.on.method${m}.init${ii}.lead${hh}.txt

m=`expr $m + 1 ` 
done
i=`expr $i + 36 `
done 
