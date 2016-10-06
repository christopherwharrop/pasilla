#!/bin/sh -l

# Load required modules
module load intel/16.1.150

# Set paths for GPTL and tools
GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

# Set the cases (method_threads) to test
cases="1_1 2_1 3_1 4_1 4_3"

# OMP settings
export OMP_NUM_THREADS=1  
export OMP_STACKSIZE=1G

# Clean out previous results
rm -f ./*.stdout*
rm -f ./*.answer*
rm -f ./timing.*

# Copy the executable
rm -f varsolver.exe
cp ../exe/varsolver.exe .

# Copy the namelist
rm -f varsolver.namelist
cp ../parm/varsolver.namelist .

# Copy the background input files
cp background/* .

# Copy the observation files
cp obs/* .

# Test each solver method
fail=0
for case in ${cases}; do

  method=`echo $case | cut -f1  -d_`
  threads=`echo $case | cut -f2  -d_`

  echo -n "Testing method ${method} using ${threads} threads... "

  # Set the number of threads for this case
  export OMP_NUM_THREADS=$threads

  # Set the method in the namelist
  sed -i "s/mthd = [[:digit:]]/mthd = ${method}/" varsolver.namelist

  # Run varsolver
  ./varsolver.exe < varsolver.namelist > varsolver.${case}.stdout

  # Post-process timing information 
  ${GPTL_PATH}/hex2name.pl ./varsolver.exe ./timing.0 > ./timing.varsolver.${case}.txt
  ${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.${case}.txt > timing.varsolver.parse.${case}.txt
  mv timing.0 timing.varsolver.${case}.hex

  # Compare the output against the baseline
  egrep 'FIN|final cost' varsolver.${case}.stdout > varsolver.${case}.answer
  basetime=`grep adept baseline/timing.varsolver.parse.${case}.txt | awk '{print $3}' `
  thistime=`grep adept timing.varsolver.parse.${case}.txt | awk '{print $3}' `
  diff varsolver.${case}.answer baseline/varsolver.${case}.answer
  if [ $? -eq 0 ]; then
    echo "PASS in ${thistime} seconds (baseline is ${basetime} seconds)"
  else
    echo "FAIL "
    fail=1
  fi

done

exit $fail
