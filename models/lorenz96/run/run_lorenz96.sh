#!/bin/sh

# Make local copies of executable and namelist
cp ../parm/lorenz96.namelist .
cp ../exe/lorenz96.exe .

# Run
./lorenz96.exe < ./lorenz96.namelist

