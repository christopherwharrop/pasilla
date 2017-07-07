#!/bin/sh -l

module load intel/16.1.150
module load netcdf
#module load pgi/16.10
#module load netcdf/4.4.0

GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp/bin
#GPTL_PATH=/contrib/gptl/gptl-v5.5_nompi_noomp_pgi/bin
PARSE_PATH=/home/Christopher.W.Harrop/bin

resol="21"

ulimit -s 2048000
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=24
export OMP_STACKSIZE=1G
#ulimit -a

# Link to the input for QG model instantiation
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/sf7910T106.shfs
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgbergT${resol}.dat
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgcoefT${resol}.dat
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qgforcingT${resol}.nc
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/qginitT${resol}.nc

# Link to the input for the QG adjoint model
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/mu_operator.dat

# Link to B matrix
ln -sf /scratch4/BMC/gsd-hpcs/QG/inputdata/bT${resol}.nc b.nc

# Link to background file
ln -sf ../test/background/QG/T${resol}/bkgin_0000002.nc

# Link to obs file
ln -sf ../test/obs/QG/qgobs_1.txt

../exe/varsolver_qg.exe < ../parm/varsolver_qg.namelist > varsolver_qg.stdout 2> varsolver_qg.stderr

${GPTL_PATH}/hex2name.pl ../exe/varsolver_qg.exe ./timing.0 > ./timing.varsolver_qg.txt

${PARSE_PATH}/parsetiming.rb -t 1 -d 6 timing.varsolver_qg.txt > timing.varsolver_qg.parse.txt
