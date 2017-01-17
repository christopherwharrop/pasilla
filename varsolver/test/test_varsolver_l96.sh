#!/bin/sh -l

# Load required modules
module purge
module load intel/16.1.150
module load netcdf/4.3.0
#module load pgi/16.10
#module load netcdf/4.4.0
module load nccmp

# Set paths for GPTL and tools
GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
#GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp_pgi/bin
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

  echo "Testing method ${method} using ${threads} threads... "

  # Set the number of threads for this case
  export OMP_NUM_THREADS=$threads

  # Set the method in the namelist
  sed -i "s/mthd = [[:digit:]]/mthd = ${method}/" ${VARSOLVER_NAMELIST}

  # Run varsolver
  ./${VARSOLVER_EXE} < ${VARSOLVER_NAMELIST} > varsolver.${case}.stdout 2> varsolver.${case}.stderr

  # Post-process timing information 
  ${GPTL_PATH}/hex2name.pl ./${VARSOLVER_EXE} ./timing.0 > ./timing.varsolver.${case}.txt
  ${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.${case}.txt > timing.varsolver.parse.${case}.txt
  mv timing.0 timing.varsolver.${case}.hex
  basetime=`grep adept ../baseline/timing.varsolver.parse.${case}.txt | awk '{print $3}' `
  thistime=`grep adept timing.varsolver.parse.${case}.txt | awk '{print $3}' `
  echo "   Completed in ${thistime} seconds (baseline is ${basetime} seconds)"

  # Compare the standard output against the baseline
  egrep 'FIN|final cost' varsolver.${case}.stdout > varsolver.${case}.answer
  cmp -s varsolver.${case}.answer ../baseline/varsolver.${case}.answer
  if [ $? -ne 0 ]; then
    fail=1
    echo "   stdout file does not match baseline"
  fi

  # Compare the NetCDF output against the baseline
#module swap pgi intel/16.1.150 2> /dev/null
#module load netcdf/4.3.0 2> /dev/null
  if [ ${method} -lt 3 ]; then
    mv bkgout_0000001.nc bkgout_0000001.${case}.nc
    nccmp -m -d -t 1e-15 bkgout_0000001.${case}.nc ../baseline/bkgout_0000001.${case}.nc
    error=$?
  else
    mv bkgout_0000002.nc bkgout_0000002.${case}.nc
    nccmp -m -d -t 1e-12 bkgout_0000002.${case}.nc ../baseline/bkgout_0000002.${case}.nc
    error=$?
  fi
  if [ $error -ne 0 ]; then
    fail=1
    echo "   NETCDF output does not match baseline"
  fi
#module swap intel pgi/16.10 2> /dev/null
#module load netcdf/4.4.0 2> /dev/null

  # Report PASS/FAIL
  if [ ${fail} -eq 1 ]; then
    echo "   FAIL"
  else
    echo "   PASS"
  fi

done
exit $fail
