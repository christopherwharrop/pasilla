#!/bin/bash -l

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

module load intel/16.1.150
module load netcdf

resol="21"

qgdir="/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.top/models/QG"
parmdir="${qgdir}/parm"
outdir="${qgdir}/outputdata"
expid='harr'
rundir="${qgdir}/rundir/run${expid}"
obsfile="sf7910T106.shfs"

plot=0

# Create an empty output directory
mkdir -p "${outdir}/$expid"
rm -f ${outdir}/$expid/*

# Create an empty run directory and cd into it
rm -rf ${rundir}
mkdir -p ${rundir}
cd ${rundir}

# Copy the inputdata into place
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/${obsfile} .
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT${resol}.dat .
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT${resol}.dat .
cp -prd /scratch4/BMC/gsd-hpcs/QG/inputdata/qgstartT${resol}.dat .

# Make sure obs file exists.
if [ ! -e $obsfile ]; then
  echo "${qgdir}/inputdata/$obsfile does not exist" ; exit 1
fi

# Copy the namelist into the run directory
cp ${parmdir}/namelist namelist.input

# Set the resolution in the namelist
sed -i "s/resolution = [[:digit:]]*/resolution = ${resol}/" namelist.input

# Set the experiment id in the namelist
sed -i "s/expid = '.*'/expid = \'${expid}\'/" namelist.input

# Set the obs file in the namelist
sed -i "s/obsfile = '.*'/obsfile = \'${obsfile}\'/" namelist.input

# Copy the executable to the run directory
cp ${qgdir}/exe/QG.exe .

# Run the model
./QG.exe

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


#rm ${ftex}.tex ${ftex}.dvi 

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

#rm ${ftex}.tex ${ftex}.dvi *.log *.aux *.eps *mean*


exit

