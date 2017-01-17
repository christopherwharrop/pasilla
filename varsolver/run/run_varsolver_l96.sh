#!/bin/sh -l

module load intel/16.1.150
#module load pgi/16.10
#module load netcdf/4.4.0

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
#GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp_pgi/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

export OMP_STACKSIZE=1G

../exe/varsolver_l96.exe < ../parm/varsolver_l96.namelist > varsolver_l96.stdout

${GPTL_PATH}/hex2name.pl ../exe/varsolver_l96.exe ./timing.0 > ./timing.varsolver_l96.txt

${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver_l96.txt > timing.varsolver_l96.parse.txt
