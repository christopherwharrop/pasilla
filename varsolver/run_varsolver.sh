#!/bin/sh -l

module load intel impi netcdf

GPTL_PATH=/contrib/gptl/gptl-v5.4.4_impi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

mpiexec ./varsolver.exe > test.out

${GPTL_PATH}/hex2name.pl ./varsolver.exe ./timing.0 > ./timing.varsolver.txt

${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.txt > timing.varsolver.parse.txt
