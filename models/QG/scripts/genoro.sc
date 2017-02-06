#!/bin/bash

# this script processes the geopotential at surface and land-sea mask at high resolution
# obtained from the ERA-INTERIM dataset at ECMWF
# the result are inputfiles with mountain height in m and land-sea mask at various
# gaussian grid resolutions

cd ../inputdata

#cdo selvar,z berg.nc z.nc
#cdo divc,9.81 z.nc height.nc
#cdo selvar,lsm berg.nc lsm.nc

#cdo remapcon2,n16 height.nc heightT21.nc
#cdo remapcon2,n32 height.nc heightT42.nc
#cdo remapcon2,n48 height.nc heightT63.nc
#cdo remapcon2,n80 height.nc heightT106.nc
#cdo remapcon2,n16 lsm.nc lsmT21.nc
#cdo remapcon2,n32 lsm.nc lsmT42.nc
#cdo remapcon2,n48 lsm.nc lsmT63.nc
#cdo remapcon2,n80 lsm.nc lsmT106.nc

for resol in T21 T42 T63 T106
do

cat > makefile.gs <<==
function main(args)
'sdfopen height${resol}.nc'
'q file'
pull return
line=sublin(result,5)
nlon=subwrd(line,3)
nlat=subwrd(line,6)

x=1
while (x <= nlon)
  'set x 'x
  y=1
  while (y <= nlat)
    'set y 'y
    'd z'
    tmp=subwrd(result,4)
    if (x=1 & y=1)
      '!echo 'tmp' > qgberg${resol}new.dat'
    else
      '!echo 'tmp' >> qgberg${resol}new.dat'
    endif
    y=y+1
  endwhile
  x=x+1
endwhile

'close 1'

'sdfopen lsm${resol}.nc'

x=1
while (x <= nlon)
  'set x 'x
  y=1
  while (y <= nlat)
    'set y 'y
    'd lsm'
    tmp=subwrd(result,4)
    '!echo 'tmp' >> qgberg${resol}new.dat'
    y=y+1
  endwhile
  x=x+1
endwhile

'quit'
return
==

grads -blc "run makefile.gs"

done





