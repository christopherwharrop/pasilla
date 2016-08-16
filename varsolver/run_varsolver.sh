#!/bin/sh -l

module load intel

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

./varsolver.exe > varsolver.stdout

${GPTL_PATH}/hex2name.pl ./varsolver.exe ./timing.0 > ./timing.varsolver.txt

${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver.txt > timing.varsolver.parse.txt
