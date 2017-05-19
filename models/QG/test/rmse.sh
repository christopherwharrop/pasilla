#!/bin/sh -l

file1=$1
file2=$2

module load intel/16.1.150
module load netcdf
module load grib_api
module load cdo

cdo=`cdo --no_warnings -s output -sqrt -fldmean -sqr -sub $file1 $file2 2>rmse.err`
if [ $? -ne 0 ]; then
  echo "ERROR: cdo FAILED!"
  cat rmse.err
  rm -f rmse.err
  exit 1
fi
echo $cdo | cut -d " " -f 3

exit 0
