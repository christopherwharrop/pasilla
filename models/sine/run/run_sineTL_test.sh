#!/bin/sh

# Make local copies of executable and namelist
cp ../parm/sine.namelist .
cp ../exe/sineTL_test.exe .

# Run
./sineTL_test.exe < ./sine.namelist

