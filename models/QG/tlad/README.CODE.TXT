CODE INFO

FROM RUNQGMODEL:

In the above, the notion is that the TL propagates the trajectory (the difference
between X(t-1) and X(t)) forward in time. That propagated trajectory, added to
the initial streamfunction (PSI), fields the streamfunction at t+1

To go the other direction, the TL propagates the trajectory forward in time one step
(in the script I have, from t=1 to t=2).  This allows for the calculation of the 
change of the trajectory. This change is then propagated back to t=0, using the AD.
That propagated change is then added onto the t=1 trajectory to get the t=0 trajectory.
After that, we add on to the t=1 state to get the t=0 state 
   
FROM THE FORTRAN CODE

The tangent linear code is in 'qgtlad-tang.F', the adjoint code is in 'qgtlad-adjo.F'.
Additionally, a batch of utility subroutines are in 'qgtlad-util.F'

Both the TL and AD are called with something like:
CALL TANG/ADJ (AS,BS,NUUR)

Where "AS" is the trajectory of input, "BS" the trajectory output, and "NUUR" something
related to the number of time steps. For now... I've only done it for a single time step.
In addition, both TANG and ADJ need the state. Instead of being passed in, this ends up being
retrieved via a subroutine "RETRIE" (after having been stored using "STOREP"). This is
not at all a good design, and will likely change. 

I know that, in the L96 code, we don't have things this way. TRAJ is an input, and then
comes out updated. Same for the STATE. As it is now, TRAJIN and TRAJOUT are two seperate 
things, and STATE is updated via the global variable route.

