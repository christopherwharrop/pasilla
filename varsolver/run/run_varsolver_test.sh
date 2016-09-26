#!/bin/sh -l

module load intel

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

export OMP_NUM_THREADS=1  
export OMP_STACKSIZE=1G
../exe/varsolver.exe < ../parm/varsolver.1.namelist > varsolver.1.stdout
${GPTL_PATH}/hex2name.pl ../exe/varsolver.exe ./timing.0 > ./timing.varsolver.1.txt
${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.1.txt > timing.varsolver.parse.1.txt
mv timing.0 timing.varsolver.1.hex

../exe/varsolver.exe < ../parm/varsolver.2.namelist > varsolver.2.stdout
${GPTL_PATH}/hex2name.pl ../exe/varsolver.exe ./timing.0 > ./timing.varsolver.2.txt
${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.2.txt > timing.varsolver.parse.2.txt
mv timing.0 timing.varsolver.2.hex

../exe/varsolver.exe < ../parm/varsolver.3.namelist > varsolver.3.stdout
${GPTL_PATH}/hex2name.pl ../exe/varsolver.exe ./timing.0 > ./timing.varsolver.3.txt
${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.3.txt > timing.varsolver.parse.3.txt
mv timing.0 timing.varsolver.3.hex

../exe/varsolver.exe < ../parm/varsolver.4.namelist > varsolver.4.stdout
${GPTL_PATH}/hex2name.pl ../exe/varsolver.exe ./timing.0 > ./timing.varsolver.4.txt
${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.4.txt > timing.varsolver.parse.4.txt
mv timing.0 timing.varsolver.4.hex

export OMP_NUM_THREADS=3
export OMP_STACKSIZE=1G

../exe/varsolver.exe < ../parm/varsolver.4.namelist > varsolver.5.stdout
${GPTL_PATH}/hex2name.pl ../exe/varsolver.exe ./timing.0 > ./timing.varsolver.5.txt
${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.5.txt > timing.varsolver.parse.5.txt
mv timing.0 timing.varsolver.5.hex
