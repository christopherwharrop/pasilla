#!/bin/sh

# Make local copies of executable and namelist
cp ../parm/sine.namelist .
cp ../exe/sineADJ_test.exe .

# Run
./sineADJ_test.exe < ./sine.namelist

