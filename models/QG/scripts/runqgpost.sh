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
module load grads

resol="21"

qgdir="/scratch4/BMC/gsd-hpcs/Christopher.W.Harrop/pasilla.top/models/QG"
parmdir="${qgdir}/parm"
outdir="${qgdir}/outputdata"
expid='harr'
rundir="${qgdir}/rundir/run${expid}"
obsfile="sf7910T106.shfs"

# cd into the run directory
cd ${rundir}

# Link to the land mask file
ln -s ${qgdir}/scripts/lpoly_mres.asc

# plot orography in file ${outdir}/$expid/qgbergT??.grads with grads
grads -blc "run ${qgdir}/scripts/oro.gs ${resol}"
rm st.gx

# plot mean streamfunction in simulation observation input with grads
grads -blc "run ${qgdir}/scripts/streamfunc_observed.gs ${resol} 800 ${obsfile}"
grads -blc "run ${qgdir}/scripts/streamfunc_observed.gs ${resol} 500 ${obsfile}"
grads -blc "run ${qgdir}/scripts/streamfunc_observed.gs ${resol} 200 ${obsfile}"
rm st.gx

# Plot mean streamfunction in simulation output with grads
grads -blc "run ${qgdir}/scripts/streamfunc_simulated.gs ${resol} 200 ${expid}"
grads -blc "run ${qgdir}/scripts/streamfunc_simulated.gs ${resol} 500 ${expid}"
grads -blc "run ${qgdir}/scripts/streamfunc_simulated.gs ${resol} 800 ${expid}"
rm st.gx

exit

