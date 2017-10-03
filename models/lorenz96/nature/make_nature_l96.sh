#!/bin/sh

# Get Lorenz96 model base directory 
L96_DIR=$(dirname $(dirname $(readlink -m $0)))

# Set runtime paramters
forcing=8.0
time_step=0.000833
one_hour=10                             # One "hour" is 10 steps
(( fcst_spinup = 40000 ))               # Spinup the model for 40,000 steps
(( fcst_length =$one_hour * 24 * 120 + $one_hour )) # Run for 120 days and one hour
output_interval=$one_hour               # Write "hourly" output

# Set up nature directory and cd into it
workdir=$L96_DIR/nature
mkdir -p $workdir
cd $workdir

# Make local copies of executable and namelist
cp ../parm/lorenz96.namelist .
cp ../exe/lorenz96.exe .

# Modify namelist parameters for spinup mode
sed -i "s/time_step = [[:digit:]\.]*/time_step = ${time_step}/" lorenz96.namelist
sed -i "s/forcing = [[:digit:]\.]*/forcing = ${forcing}/" lorenz96.namelist
sed -i "s/start_step = [[:digit:]]*/start_step = 0/" lorenz96.namelist
sed -i "s/run_steps = [[:digit:]]*/run_steps = ${fcst_spinup}/" lorenz96.namelist
sed -i "s/output_interval_steps = [[:digit:]]*/output_interval_steps = ${fcst_spinup}/" lorenz96.namelist

# Spinup the model
echo "Spinning up the model for ${fcst_spinup} steps"
./lorenz96.exe < ./lorenz96.namelist

# Change the startoutput time for run mode
sed -i "s/start_step = [[:digit:]]*/start_step = ${fcst_spinup}/" lorenz96.namelist
sed -i "s/run_steps = [[:digit:]]*/run_steps = ${fcst_length}/" lorenz96.namelist
sed -i "s/output_interval_steps = [[:digit:]]*/output_interval_steps = ${output_interval}/" lorenz96.namelist

# Run the model
echo "Running the model for ${fcst_length} steps"
./lorenz96.exe < ./lorenz96.namelist
