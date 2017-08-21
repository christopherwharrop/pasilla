Verification of 3-level QG Model

The process is, from the present directory, is to run this:
verify_qg_forecasts.sh

This makes 11 different directories, one for each lead time from 000 to 360.
It then copies in 'run_verify.sh' into each directory.

It then runs for these 11 lead times in each directory, simultaneaously.
The parallel runs are needed because it takes so long to run.

One bit: verify.exe is needed, and that comes from a source file in varsolver/src/verify_qg.f90
