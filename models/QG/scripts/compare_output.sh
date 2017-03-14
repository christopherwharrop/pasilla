#!/bin/sh

dir1=$1
dir2=$2

for ncfile in `(cd $dir1; ls -1 *.nc; cd ..)`; do
  if [ -f $dir2/$ncfile ]; then
    nccmp -d -m $dir1/$ncfile $dir2/$ncfile
    if [ $? -ne 0 ]; then
      exit
    fi
  fi
done

for file in `(cd $dir1; ls -1 | grep -v '\.nc' | grep -v '\.exe' | grep -v timing; cd ..)`; do
  if [ -f $dir2/$file ]; then
    diff -qs $dir1/$file $dir2/$file
    if [ $? -ne 0 ]; then
      exit
    fi
  fi
done
