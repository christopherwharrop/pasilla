#!/bin/sh

# Make local copies of executable and namelist
cp ../parm/lorenz96.namelist .
cp ../exe/lorenz96TL_test.exe .

# Run
./lorenz96TL_test.exe < ./lorenz96.namelist

