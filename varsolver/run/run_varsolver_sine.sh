#!/bin/sh -l

module load intel/16.1.150

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

export OMP_STACKSIZE=1G

../exe/varsolver_sine.exe < ../parm/varsolver_sine.namelist > varsolver_sine.stdout

${GPTL_PATH}/hex2name.pl ../exe/varsolver_sine.exe ./timing.0 > ./timing.varsolver_sine.txt

${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver_sine.txt > timing.varsolver_sine.parse.txt
