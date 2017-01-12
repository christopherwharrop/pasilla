#!/bin/sh -l

# Load required modules
module load intel/16.1.150
module load netcdf
module load nccmp

# Set paths for GPTL and tools
GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

# Set the cases (method_threads) to test
cases="1_1 2_1 3_1 4_1 4_3"

# OMP settings
export OMP_NUM_THREADS=1  
export OMP_STACKSIZE=1G

# Clean out previous results
workdir=./work
rm -rf $workdir
mkdir -p $workdir
cd $workdir

# Link the executable
VARSOLVER_EXE=varsolver_l96.exe
rm -f ${VARSOLVER_EXE}
ln -s ../../exe/${VARSOLVER_EXE}

# Copy the namelist
VARSOLVER_NAMELIST=varsolver_l96.namelist
rm -f ${VARSOLVER_NAMELIST}
cp ../../parm/${VARSOLVER_NAMELIST} .

# Copy the background input files
cp -prd ../background/bkgin* .
cp -prd ../background/lorenz96* .

# Copy the observation files
cp ../obs/lorenz* .

# Test each solver method
fail=0
for case in ${cases}; do

  method=`echo $case | cut -f1  -d_`
  threads=`echo $case | cut -f2  -d_`

  echo "Creating new baseline for method ${method} using ${threads} threads... "

  # Set the number of threads for this case
  export OMP_NUM_THREADS=$threads

  # Set the method in the namelist
  sed -i "s/mthd = [[:digit:]]/mthd = ${method}/" ${VARSOLVER_NAMELIST}

  # Run varsolver
  ./${VARSOLVER_EXE} < ${VARSOLVER_NAMELIST} > varsolver.${case}.stdout

  # Post-process timing information 
  ${GPTL_PATH}/hex2name.pl ./${VARSOLVER_EXE} ./timing.0 > ./timing.varsolver.${case}.txt
  ${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.${case}.txt > ../baseline/timing.varsolver.parse.${case}.txt
  mv timing.0 timing.varsolver.${case}.hex

  # Compare the standard output against the baseline
  egrep 'FIN|final cost' varsolver.${case}.stdout > ../baseline/varsolver.${case}.answer

  # Compare the NetCDF output against the baseline
  if [ ${method} -lt 3 ]; then
    mv bkgout_0000001.nc ../baseline/bkgout_0000001.${case}.nc
  else
    mv bkgout_0000002.nc ../baseline/bkgout_0000002.${case}.nc
  fi

done

exit