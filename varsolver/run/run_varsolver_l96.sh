#!/bin/sh -l

module load intel/16.1.150

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

export OMP_STACKSIZE=3G

../exe/varsolver_l96.exe < ../parm/varsolver.namelist > varsolver.stdout

${GPTL_PATH}/hex2name.pl ../exe/varsolver_l96.exe ./timing.0 > ./timing.varsolver.txt

${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.txt > timing.varsolver.parse.txt
