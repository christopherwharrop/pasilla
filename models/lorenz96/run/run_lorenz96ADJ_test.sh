#!/bin/sh

# Make local copies of executable and namelist
cp ../parm/lorenz96.namelist .
cp ../exe/lorenz96ADJ_test.exe .

# Run
./lorenz96ADJ_test.exe < ./lorenz96.namelist

