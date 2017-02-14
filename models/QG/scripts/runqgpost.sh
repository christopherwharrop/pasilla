#!/bin/bash -l

# This script runs the GrADS graphics scripts.

module load grads

# Set the resolution and experiment id
resol="21"
expid='harrop'

# Set the path to the QG model install - DO NOT CHANGE THIS LINE
qgdir=$(dirname $(dirname $(readlink -m $0) ))

# Set the path to the run directory
rundir="${qgdir}/run/${expid}"

# Set the name of the obs file to use - This must always be the T106 file (for now).
obsfile="sf7910T106.shfs"

# cd into the run directory
cd ${rundir}

# Link to the land mask file
ln -fs ${qgdir}/scripts/lpoly_mres.asc

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

