UPDATED TANGENT LINEAR AND ADJOINT CODE WORK

****
runqgmodel.sc runs the QG model forward 8 time steps, and also runs the TL and AD
for each of those time steps. The output ends up in 'outputdir*0000', with *=T
being the true TL and AD, *=X being 'fake'.

Output is in the file qgmodelT106.grads, accessed via grads and qgmodelT106.ctl
Output, frames 1-25 of that file (times are 'fake') are:

TIME =  1  2  3  4  5  6  7  8
-----  -- -- -- -- -- -- -- -- 
IC =    1  4  7 10 13 16 19 22
TL =    2  5  8 11 14 17 20 23
AD =    3  6  9 12 15 18 21 24
FW =    4  7 10 13 16 19 22 25

The TL is the TL started from the IC. The AD is the AD started from the output of the TL.
The FW is the true model run forward one time step, started from IC.

If perfect, the TL and FW frames for each of the 8 times would be identical.
If perfect, the AD and IC frames for each of the 8 times would be identical.

Output of log files and such is stored in 'rundir'

****
runqgmodel-full.sc does the same thing, but offers more detailed output. Though only
for 3 time steps, rather than 8, the output in outputdir/*0001 consists of:

TIME =  	 1    2    3
----------	--   --   --
X(t) = IC	 1    9   17 
X'(t)		 2   10   18   
X'(t+1))	 3   11   19
X(t+1) via TL    4   12   20
X'(t+2)	         5   13   21
X"(t+1)          6   14   22
X'(t)	         7   15   23
X(t) via AD	 8   16   24
X(t+1) via FW    9   17   25

As with the first script, if perfect...
The AD and IC fields will be the same
The TL and FW fields will be the same

But in addition
The two X(t) fields, '2 and 7', ought to be the same
The rest of the fields help diagnose why the above is or is not the case
Running the above scripts will produce the desired fields.
