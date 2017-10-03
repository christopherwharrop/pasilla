#!/bin/sh

# Get Lorenz96 model base directory 
L96_DIR=$(dirname $(dirname $(readlink -m $0)))

# Make obs for methods 1,2,3
for m in `seq 1 3`; do

  echo "Making obs for method $m"

  # Create the obs for method m
  ./create_obs.py $m

  # Set up obs directory and move obs into it
  obsdir=method_$m
  rm -rf $obsdir
  mkdir -p $obsdir
  mv lorenz96*.txt $obsdir

done

# Make links for methods 4 and 5
rm -f method_4 method_5
ln -s method_3 method_4
ln -s method_3 method_5
