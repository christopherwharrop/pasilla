#!/bin/bash

# this script produces a main fortran program that links with the ../src/qgmodel.F routines
# to calculate the forcing from a datafile with observations (given by obsfile)
# and to integrate the model at a resolution of T21, T42, T63 or T106 (set by resol)
# the outputdata is stored in directory ${outdir}/$expid
# the model is configured with a namelist file that is produced in this script
# the outputdata is plotted using the grads package (grads.iges.org)
# latex is used to combine the plots in a single pdf file

# resol sets the resolution
# qgdir must point to the directory of the distribution
# expid identifies the directory where the output of the run will be stored
# obsfile is the datafile with observations used to calculate the forcing from

resol="106"
GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin
compiler=ifort
GPTLFLAGS='-I/contrib/gptl/gptl-v5.5_nompi_noomp/include -L/contrib/gptl/gptl-v5.5_nompi_noomp/lib -lgptl'
fflags='-extend_source 132 -g -traceback -O2 -mcmodel medium -shared-intel -convert big_endian -finstrument-functions'
## fflags='-extend_source 132 -O0 -g -traceback -check all -fpe0 -mcmodel medium -shared-intel -convert big_endian -finstrument-functions'
qgdir="/scratch3/BMC/gsd-hpcs/Brian.Etherton/superQG-TLAD/"
outdir="${qgdir}/outputdata"
expid='0000'
rundir="${qgdir}/rundir/run${expid}"
obsfile="sf7910T106.shfs"

plot=0

if [ ${plot} == 0 ]; then


if [ -e ${outdir}/$expid ]; then
  echo "${outdir}/$expid already exists" 
  read -p "Are you sure to continue ? [YN]" -n 1 -r
  echo 
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  fi
fi

if [ ! -e ${qgdir}/inputdata/$obsfile ]; then
  echo "${qgdir}/inputdata/$obsfile does not exist" ; exit 1
fi

mkdir -p "${outdir}/$expid"
rm -f ${outdir}/$expid/*
mkdir -p ${rundir}
cd ${rundir}

# now produce the namelist file to configure the model
#
# expid = four digit number used to specify the outputdata directory
# inf = logical if true then the forcing is read from file inputdata/qgpvforT??.dat
# obsf = logical if true then the forcing is calculated from obsfile
# readstart = logical if true initial state is read from qgstartT??.dat else state of rest
# nstepsperday = number of timesteps the model takes for one day
# nstepsbetweenoutput = number of timesteps between the outputted states to the datafile
# ndayskip = length of the transient integration in days; model starts from inputdata/qgstartT??.dat
# nday = length of the trajectory in days that is outputted to the datafile
# obsfile = inputdata/$obsfile contains the observed states from which to calculate the forcing
#
# tdis= time scale in days of the Ekman damping (linear damping on vorticity at lowest level)
# addisl= parameter of the land-seamask dependent Ekman damping (more friction over land)
# addish= parameter of the orography dependent Ekman damping (more friction over slopes)
# trel= time scale of the temperature cooling in days
# tdif= time scale in days of the scale selective damping at the three levels for the smallest wavenumber
# idif= power of the laplacian for the scale selective damping 
# h0= scale height of the topography
# rrdef1= Rossby radius of deformation of the 200-500 hpa layer
# rrdef2= Rossby radius of deformation of the 500-800 hPa layer
#
# NOTE: depending on resolution the scale selective damping should be set differently as
# well as the nstepsperday
#
# 


cat > namelist <<==
 &control
 expid = '$expid',
 inf = .false.,
 obsf = .true.,
 readstart=.true.,
 nstepsperday = 72,
 nstepsbetweenoutput = 3,
 ndayskip = 10,
 nday = 1,
 obsfile = '$obsfile' ,
 /
 &param
 tdis=5.0,
 addisl=0.6,
 addish=0.4,
 trel=25.,
 tdif=3,
 idif=2,
 h0=6.,
 rrdef1=.110,
 rrdef2=.070,
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
c     nsh23:  nsh2*nvl (assumed to be 3)
c     ntadj:  number of time-steps for the TL and AD
c     ngp:    number of grid points of the Gaussian grid
c     nmat:   for the adjoint
c 
      integer nm,nlon,nlat,nvl,ntl,nsh,nsh2,nsh23,ntadj,ngp,nmat
==

# https://climatedataguide.ucar.edu/climate-model-evaluation/common-spectral-model-grid-resolutions
# if [ $resol == 85 ]; then
# cat >> truncation.h <<==
#       character*2 ft
#       parameter ( nm=85, nlon=256, nlat=128, nvl=3, ntl=nvl-1, ft="85")
#       parameter ( nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)
# ==
# fi

if [ $resol == 106 ]; then
cat >> truncation.h <<==
      character*3 ft
      parameter ( nm=106,nlon=320,nlat=160,nvl=3,ntl=nvl-1,ft="106")
      parameter ( nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)
      parameter ( nsh23=NSH2*NVL, NTADJ=12 )
      parameter ( NMAT= 22*23-22) 
==
fi

if [ $resol == 63 ]; then
cat >> truncation.h <<==
      character*2 ft
      parameter ( nm=63, nlon=192, nlat=96, nvl=3, ntl=nvl-1, ft="63")
      parameter ( nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)
      parameter ( nsh23=NSH2*NVL, NTADJ=12 )
      parameter ( NMAT= 22*23-22)
==
fi

if [ $resol == 42 ]; then
cat >> truncation.h <<==
      character*2 ft
      parameter ( nm=42, nlon=128, nlat=64, nvl=3, ntl=nvl-1, ft="42")
      parameter ( nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)
      parameter ( nsh23=NSH2*NVL, NTADJ=12 )
      parameter ( NMAT= 22*23-22)
==
fi

if [ $resol == 21 ]; then
cat >> truncation.h <<==
      character*2 ft
      parameter ( nm=21, nlon=64, nlat=32, nvl=3, ntl=nvl-1, ft="21")
      parameter ( nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh, ngp=nlon*nlat)
      parameter ( nsh23=NSH2*NVL, NTADJ=12 )
      parameter ( NMAT= 22*23-22)
==
fi

cp ${qgdir}/src/qgtlad*F .
cp ${qgdir}/src/comqg*h .
cp ${qgdir}/src/qgmodel*F .
cp ${qgdir}/src/nag.f .

cat > runqgmodel.F <<==
c23456789012345678901234567890123456789012345678901234567890123456789012
      program runqgmodel
c-----------------------------------------------------------------------
c *** integrates the qgmodel with parameters from inputdata/namelist
c-----------------------------------------------------------------------
      use gptl
      implicit none
     
#include "truncation.h"
#include "comqg.h"
      
      integer istep,itlad,nstep 
      integer ret
      real*8  oldq(nsh2,nvl),olds(nsh2,nvl),oldt(nsh2,ntl)
      real*8  newq(nsh2,nvl),news(nsh2,nvl),newt(nsh2,ntl)
      real*8  tngt(nsh2,nvl),tpls(nsh2,nvl)
      real*8  tlin(nsh2,nvl),tlout(nsh2,nvl),ys(nsh2,nvl)

      ret = gptlinitialize ()                     ! Initialize GPTL
      ret = gptlstart ('qgmodel')                 ! Start a timer     
      ret = gptlstop ('qgmodel')                  ! stop timer for the entire code execution
 
      rootdir="${qgdir}"
			
      call initqg
      write(*,*) 'Experiment ',expid
      write(*,*) 
      write(*,*) 'Integrating transient days: ',ndayskip
      
      nstep=ndayskip/dt
			
      do istep=1,nstep
        call forward
	CALL QTOPSI
        if(istep.gt.1) then
           oldq = newq
           olds = news
           oldt = newt
        else
           oldq = 0.0
           olds = 0.0
           oldt = 0.0
        endif
	news = psi
	newt = psit
	newq = qprime
      enddo
			
      write(*,*) 'Integrating trajectory of days: ',nday
      
      istep=0
      nstep=nday/dt
      call diag(0)

44    FORMAT(4F14.8)
c     EDIT FROM BJE
      do istep=1,8
	print *,"ISTEP",istep

C       STORE THE STATE OF THE SYSTEM FOR TL/AD, AND FUTURE RUN OF FULL MODEL
	psi    = news
        call storep(0)

C	FIRST, RUN TL FROM QPRIME STATE
  	dtt    = dt*pi*4d0
        tlin   = news-olds
	tngt   = tlin

C       USE TL TO ADVANCE THE TRAJ 1 STEP
        psi    = news
        call psitoq
        call tang(tlin,tlout,0)
	tngt   = tlout

C       PLOT OF THE NEW STATE (news+traj)
C	X (t+1) = X (t) + X'(t+1)
        tpls   = news+tngt
        psi    = tpls
	call psitoq
        call diag(nstepsbetweenoutput)

C       SECOND, RUN ADJOINT FROM TL STATE - SHOULD END UP BACK AT THE START 
C       COMPARE THIS FIELD TO QPRIME FROM PRIOR TIME

C       PROPAGATE THE TRAJECTORY BACKWARD ONE STEP
	call storep(0)
	call psitoq
	dtt    = -dt*pi*4d0
 	tlin   = tngt
	call tang(tlin,tlout,0)
	tngt   = tlout

	tlin = tlin-tlout
	psi  = news
	dtt    = dt*pi*4d0
	call adj(tlin,tlout,0)

C       PUT TOGETHER THE STREAMFUNCTION FIELD
C	THIS IS X (t) = X (t+1) - X'(t+1)
C		X'(t) = X'(t+1) - X"(t+1)
	tngt   = tngt-tlout
  	psi    = tpls-tngt-2.*tlout
        call psitoq
        call diag(nstepsbetweenoutput)

C       THIRD, RUN FORWARD MODEL FROM QPRIME - GET NEXT TIME STEP
C	COMPARE THIS FIELD TO THE FIRST
 	dtt    = dt*pi*4d0
	psi    = news
	psit   = newt
	qprime = newq
	print *,"IN FORWARD, DTT=",dtt
        call forward
        call diag(nstepsbetweenoutput)

C       FOURTH, UPDATE FIELDS FOR NEXT GO
	oldt   = newt
	olds   = news
        oldq   = newq
        newt   = psit
        news   = psi
        newq   = qprime 

      enddo
			
      call writestate

      ret = gptlpr (0)                         ! print timers for MPI rank 0
      ret = gptlfinalize ()
 
c     return
      end
==

$compiler $fflags -c  qgmodel.F -o qgmodel.o
$compiler $fflags -c  nag.f -o nag.o

$compiler $fflags -c qgmodel-tlad.F -o qgmodel-tlad.o
$compiler $fflags -c qgtlad-util.F -o qgtlad-util.o
$compiler $fflags -c qgtlad-tang.F -o qgtlad-tang.o
$compiler $fflags -c qgtlad-adjo.F -o qgtlad-adjo.o

$compiler $fflags -I${qgdir}/src -o runqgmodel runqgmodel.F qgmodel-tlad.o qgmodel.o nag.o $GPTLFLAGS
$compiler $fflags -I${qgdir}/src -o runqgtlad  runqgmodel.F qgtlad-adjo.o qgtlad-tang.o  qgtlad-util.o qgmodel.o nag.o $GPTLFLAGS

echo "RUNQGMODEL"    
./runqgmodel > runqgmodel.log 
${GPTL_PATH}/hex2name.pl ./runqgmodel ./timing.0 > ./timing.qgorig.txt
${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.qgorig.txt > timing.qgorig.parse.txt
mv ${outdir}/$expid ${outdir}/X$expid
mkdir ${outdir}/$expid

echo "RUNQGTLAD"
./runqgtlad > runqgtlad.log
${GPTL_PATH}/hex2name.pl ./runqgtlad ./timing.0 > ./timing.qgtlad.txt
${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.qgtlad.txt > timing.qgtlad.parse.txt

# rm *.o runqgmodel.F runqgmodel
fi


cd ${outdir}/$expid

# plot orography in file ${outdir}/$expid/qgbergT??.grads  with grads

cat > plot.gs <<==1

function test(args)

offset=0
* light grey (reserve 81-89)
'set rgb 81 180 180 180'
'set rgb 82 160 160 160'

'open qgbergT'${resol}'.ctl'

'set lon -180 180'
'set lat 0 90'
'set vpage 0 11 0 8.5'
'set mproj nps'

'set grads off'
'set font 0'
'set gxout contour'

'set xlopts 1 3 0.18'
'set ylopts 1 3 0.18'

'set clevs 9e32'
'set mpdraw off'
'set grid off'
'd oro'
'run basemap.gs L 81 81 M'

'set ccolor 1'
'set cthick 7'
'set clopts 1 1.2 0.1'
'set clab off'
'set mpdraw off'
'set grid on 5 82 10'
'set poli off'
'set xaxis 0 360 30'
'set ylevs 15 30 45 60 75'
'd oro'
'set string 1 tc'
'set strsiz 0.16 0.18'
'draw string 5.5 8. ECMWF orography at T'$resol

'enable print st.gx'
'print'
'disable print'
'!gxeps -c st.gx'
'!epstopdf st.eps'
'!mv st.pdf oro${resol}.pdf'
'!mv st.eps oro${resol}.eps'
'quit'
return
==1

grads -blc "run plot.gs"

rm st.gx

# plot mean streamfunction in file inputdata/obsfile  with grads

cat > plot.gs <<==1
function test(args)

lev=subwrd(args,1)

* light grey (reserve 81-89)
'set rgb 81 180 180 180'
'set rgb 82 160 160 160'

'open ${qgdir}/inputdata/${obsfile}${resol}.ctl'

pi=3.14159265
radius=6.37e+6
'd 4.0*'pi'/(24.0*3600.0)'
om=subwrd(result,4)
'd 'radius'*'radius'*'om'*1e-7'
scalesf=subwrd(result,4)

'q file'
line=sublin(result,5)
ntime=subwrd(line,12)

'set lon -180 180'
'set lat 0 90'
'set lev 'lev
'set vpage 0 11 0 8.5'
'set mproj nps'

'set grads off'
'set font 0'
'set gxout contour'

'set xlopts 1 3 0.18'
'set ylopts 1 3 0.18'

'set clevs 9e32'
'set mpdraw off'
'set grid off'
'd sf'
'run basemap.gs L 81 81 M'

'set ccolor 1'
'set cthick 7'
'set clopts 1 1.2 0.1'
*'set clab off'
'set mpdraw off'
'set grid on 5 82 10'
'set poli off'
'set xaxis 0 360 30'
'set ylevs 15 30 45 60 75'
'd 'scalesf'*ave(sf,t=1,t='ntime')'
'set string 1 tc'
'set strsiz 0.16 0.18'
'draw string 5.5 8. Mean ${obsfile} streamfunction at 'lev

'enable print st.gx'
'print'
'disable print'
'!gxeps -c st.gx'
'!epstopdf st.eps'
'!mv st.pdf ${obsfile}.mean'lev'.pdf'
'!mv st.eps ${obsfile}.mean'lev'.eps'
'quit'
return
==1

grads -blc "run plot.gs 800"
grads -blc "run plot.gs 500"
grads -blc "run plot.gs 200"

rm plot.gs st.gx

# combine the plots in a single pdf file with latex

ftex="${obsfile}"

cat > ${obsfile}.tex <<==
\documentclass[]{article}
\usepackage{graphicx}
\usepackage{mathptmx}
\pagestyle{empty}
\setlength{\textheight}{234mm}
\setlength{\textwidth}{17.4cm}
\setlength{\oddsidemargin}{-1cm}
\begin{document}
\begin{figure}
\begin{center}
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./${obsfile}.mean200.eps}
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./${obsfile}.mean500.eps}\\
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./${obsfile}.mean800.eps}
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./oro${resol}.eps}
\end{center}
\end{figure}
\end{document}
==
# latex ${ftex}.tex
# latex ${ftex}.tex
# dvipdf ${ftex}.dvi


rm ${ftex}.tex ${ftex}.dvi 

cd ${outdir}/$expid

# plot mean streamfunction of the integration in the datafile qgmodelsfT??.grads

cat > plot.gs <<==1

function test(args)

lev=subwrd(args,1)
year=subwrd(args,2)

offset=0
* light grey (reserve 81-89)
'set rgb 81 180 180 180'
'set rgb 82 160 160 160'

'open qgmodelsfT${resol}.ctl'

'q file'
line=sublin(result,5)
ntime=subwrd(line,12)

'set lon -180 180'
'set lat 0 90'
'set lev 'lev
'set vpage 0 11 0 8.5'
'set mproj nps'

'set grads off'
'set font 0'
'set gxout contour'

'set xlopts 1 3 0.18'
'set ylopts 1 3 0.18'

'set clevs 9e32'
'set mpdraw off'
'set grid off'
'd psi'
'run basemap.gs L 81 81 M'

'set ccolor 1'
'set cthick 7'
'set clopts 1 1.2 0.1'
*'set clab off'
'set mpdraw off'
'set grid on 5 82 10'
'set poli off'
'set xaxis 0 360 30'
'set ylevs 15 30 45 60 75'
'd 1e-7*ave(psi,t=1,t='ntime')'
'set string 1 tc'
'set strsiz 0.16 0.18'
'draw string 5.5 8. Mean T${resol}-${expid} streamfunction at 'lev

'enable print st.gx'
'print'
'disable print'
'!gxeps -c st.gx'
'!epstopdf st.eps'
'!mv st.pdf T${resol}-${expid}.mean'lev'.pdf'
'!mv st.eps T${resol}-${expid}.mean'lev'.eps'
'quit'
return
==1

grads -blc "run plot.gs 800 "
grads -blc "run plot.gs 500 "
grads -blc "run plot.gs 200 "

rm plot.gs st.gx

cd ${outdir}/$expid

# combine the observed and simulated mean streamfunctions in one pdf file with latex

ftex="T${resol}-${expid}"

cat > ${ftex}.tex <<==
\documentclass[]{article}
\usepackage{graphicx}
\usepackage{mathptmx}
\pagestyle{empty}
\setlength{\textheight}{234mm}
\setlength{\textwidth}{17.4cm}
\setlength{\oddsidemargin}{-1cm}
\begin{document}
\begin{figure}
\begin{center}
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./${obsfile}.mean200.eps}
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./T${resol}-${expid}.mean200.eps}\\
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./${obsfile}.mean500.eps}
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./T${resol}-${expid}.mean500.eps}\\
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./${obsfile}.mean800.eps}
\includegraphics[width=0.45\textwidth,angle=-90,clip]
{./T${resol}-${expid}.mean800.eps}
\end{center}
\end{figure}
\end{document}
==
# latex ${ftex}.tex
# latex ${ftex}.tex
# dvipdf ${ftex}.dvi

rm ${ftex}.tex ${ftex}.dvi *.log *.aux *.eps *mean*

mv ${outdir}/$expid ${outdir}/T$expid

exit

