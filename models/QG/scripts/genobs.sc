#!/bin/bash
#
# this script processes the netcdf file downloaded from ERA-INTERIM relvor7910.nc
# which contains six hourly relative vorticity data from DJF 1979-2010 at
# 800, 500 and 200 hPa
# the result is a file with streamfunction spherical harmonical coefficients at
# T106 resolution that is read by the qgmodel to determine the forcing
# the name of that file is sf7910T106.shfs

cd ../inputdata
#cdo remapbil,n80 relvor7910.nc relvor7910n80.nc
#cdo daymean relvor7910n80.nc relvor7910dn80.nc

cat > relvorn80.ctl <<==
DSET ^relvor7910dn80.nc
TDEF time 11192 LINEAR 0Z01DEC1979 1DY
==


cat > makefile.gs <<==
function main(args)
'xdfopen relvorn80.ctl'
'q file'

line=sublin(result,5)
nlon=subwrd(line,3)
nlat=subwrd(line,6)

'set x 1 'nlon
'set y 1 'nlat

'set gxout fwrite'
'set fwrite -be -sq  relvor7910n80.grads'

t=1

while (t <= 2798)
  'set t 't
  'set lev '800
  'd vo'
  'set lev '500
  'd vo'
  'set lev '200
  'd vo'
  t=t+1
endwhile

'quit'
return
==

grads -blc "run makefile.gs"

rm makefile.gs

resol="106"
compiler='gfortran'
fflags="-fbounds-check -O2 -ffixed-line-length-none -fconvert=big-endian"
qgdir="/Users/seltini/Werk/qgmodel42"
rundir="${qgdir}/rundir/genobs"

mkdir -p ${rundir}
cd ${rundir}

cat > namelist <<==
 &control
  inf = .false.,
/
 &param
 /
==
cat > truncation.h <<==
c *** PARAMETERS
c     nm  :   the truncation is of type T(riangular) nm. 
c     nlon:   number of longitude points of the Gaussian grid
c     nlat:   number of latitude  points of the Gaussian grid
c     nvl :   number of vorticity levels in the vertical 
c             (should be set to 3)
c     ntl :   number of temperature levels in the vertical 
c             (equal to nvl-1)
c     nsh :   half of nsh2
c     nsh2:   number of coefficients needed to define one level of the 
c             T nm model
c     ngp:    number of grid points of the Gaussian grid
c 
      integer nm,nlon,nlat,nvl,ntl,nsh,nsh2,ngp
==

if [ $resol == 106 ]; then
cat >> truncation.h <<==
      character*3 ft
      parameter ( nm=106, nlon=320, nlat=160, nvl=3, ntl=nvl-1, ft="106")
      parameter ( nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)
==
fi

cp ${qgdir}/src/comqg.h .
cp ${qgdir}/src/qgmodel.F .
cp ${qgdir}/src/nag.f .

cat > genobs.f <<==
c23456789012345678901234567890123456789012345678901234567890123456789012
      program genobs
c-----------------------------------------------------------------------
c *** calculates non-dimensional streamfunction observation file
c-----------------------------------------------------------------------
      implicit none
      
#include "truncation.h"
#include "comqg.h"
      
      integer istep,nstep,i,j,l,nvar
      real*4  vo4(nlat,nlon)
      real*8  vosp(nsh2),sfsp(nsh2,nvl),vo(nlat,nlon),facnondim
      real*8  sfspfs(nsh2,nvl)

      
      nvar=(nm+2)*nm
      
      rootdir="${qgdir}"
			
      call initqg
      facnondim=24d0*3600d0/(4d0*pi)
      
      open(unit=21,
     *file=rootdir(1:rootdirl)//'/inputdata/relvor7910n80.grads',
     *form='unformatted')
      open(unit=22,
     *file=rootdir(1:rootdirl)//'/inputdata/sf7910T106.shfs',
     *form='unformatted')
      
      nstep=2798
      
      do istep=1,nstep
        do l=1,nvl
          read(21)  ((vo4(j,i),i=1,nlon),j=1,nlat)
          do i=1,nlon
            do j=1,nlat
              vo(j,i)=vo4(j,i)*facnondim
            enddo
          enddo
          call ggtosp(vo,vosp)
          call lapinv(vosp,sfsp(1,nvl-l+1))
        enddo        
        call fmtofs(sfsp,sfspfs)
        do l=nvl,1,-1
          write(22) (real(sfspfs(i,l)),i=1,nvar)
        enddo
      enddo 
      
      return
      end
==

$compiler $fflags -c  qgmodel.F -o qgmodel.o
$compiler $fflags -c  nag.f -o nag.o

$compiler $fflags -I${qgdir}/src -o genobs genobs.F qgmodel.o nag.o

./genobs 

exit