#!/bin/sh

# TESTING TO ENSURE THAT THE VARSOLVER RUNS ARE AS THEY SHOULD BE
./run_varsolver.sh
# OMP_STACKSIZE=1G
# CASE 1: NUM_OMP_THREADS=1:  3DVAR, OBS AT TIMES 1 AND 3 COMPARED TO BACKROUND AT TIME 2
# CASE 2: NUM_OMP_THREADS=1:  3DVAR, OBS MATCHED TO BACKGROUND CORRECTLY
# CASE 3: NUM_OMP_THREADS=1:  SEQUENTIAL 4DVAR
# CASE 4: NUM_OMP_THREADS=1:  TIME-PARALLEL 4DVAR
# CASE 5: NUM_OMP_THREADS=3:  TIME-PARALLEL 4DVAR

# TEST FOR COST FUNCTION CORRECTNESS - SHOULD BE SAME TO 6 DECIMAL PLACES OR SO
grep 'final cost' varsolver.*.stdout
# varsolver.1.stdout: final cost =    72.6101726236748       after           51  iterations
# varsolver.2.stdout: final cost =    49.3304146396658       after           13  iterations
# varsolver.3.stdout: final cost =    48.9378665525746       after           13  iterations
# varsolver.4.stdout: final cost =    48.9414864060254       after           13  iterations
# varsolver.5.stdout: final cost =    48.9414864060254       after           13  iterations

# TEST FOR FINAL VALUES - SHOULD BE THE SAME TO 3 DECIMAL PLACES OR SO - TWO GRIDPOINTS
grep "FIN    2  965" varsolver.*.stdout
# varsolver.1.stdout:     FIN    2  965   55.1758   55.4867   52.3553
# varsolver.2.stdout:     FIN    2  965   52.3587   55.4867   52.3553
# varsolver.3.stdout:     FIN    2  965   52.4434   55.4867   52.3553
# varsolver.4.stdout:     FIN    2  965   52.3989   55.4867   52.3553
# varsolver.5.stdout:     FIN    2  965   52.3989   55.4867   52.3553
grep "FIN    2 3000" varsolver.*.stdout
# varsolver.1.stdout:     FIN    2 3000   44.6490   49.9999   46.8605
# varsolver.2.stdout:     FIN    2 3000   46.9991   49.9999   46.8605
# varsolver.3.stdout:     FIN    2 3000   46.9517   49.9999   46.8605
# varsolver.4.stdout:     FIN    2 3000   47.0007   49.9999   46.8605
# varsolver.5.stdout:     FIN    2 3000   47.0007   49.9999   46.8605

# TEST FOR TIMING VALUES - SHOULD NOT VARY TOO MUCH
grep 'adept' timing.varsolver.parse.*.txt
# timing.varsolver.parse.1.txt:  adept                                      1    127.6    100.0%
# timing.varsolver.parse.2.txt:  adept                                      1     35.0    100.0%
# timing.varsolver.parse.3.txt:  adept                                      1    108.8    100.0%
# timing.varsolver.parse.4.txt:  adept                                      1    108.6    100.0%
# timing.varsolver.parse.5.txt:  adept                                      1     44.2    100.0% 
