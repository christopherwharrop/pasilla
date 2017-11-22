#!/bin/sh

# Get Lorenz96 model base directory 
QG_DIR=$(dirname $(dirname $(readlink -m $0)))

# Create the obs for method 1,2,3
./create_obs.py

 # Set up obs directory and move obs into it
for m in `seq 1 3`; do
  obsdir=method_$m
  rm -rf $obsdir
  mkdir -p $obsdir
  mv qgobs_${m}_*.txt $obsdir
done

# Make links for methods 4 and 5
rm -f method_4 method_5
ln -s method_3 method_4
ln -s method_3 method_5
