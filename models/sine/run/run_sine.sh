#!/bin/sh

# Make local copies of executable and namelist
cp ../parm/sine.namelist .
cp ../exe/sine.exe .

# Run
./sine.exe < ./sine.namelist

