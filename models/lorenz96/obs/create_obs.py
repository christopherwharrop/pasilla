#!/bin/env python

import sys
from scipy.io import netcdf
import numpy as np

method=int(sys.argv[1])
print "Generating obss for method = ", method

text=open("obs_sites.txt",'r')
lobs=text.readlines()
print len(lobs)
text.close()

spinup = 40000
assim_window=10

#for h in range(spinup + 120,spinup + 28801,120):
for h in range(spinup + 60,spinup + 28801,60):
    fnm_tm1 ="../nature/lorenz96out_"+str(h-assim_window).zfill(7)+".nc"
    fnm_t   ="../nature/lorenz96out_"+str(h).zfill(7)+".nc"
    fnm_tp1 ="../nature/lorenz96out_"+str(h+assim_window).zfill(7)+".nc"

    print h,fnm_tm1, fnm_t, fnm_tp1

    ncf_tm1 = netcdf.netcdf_file(fnm_tm1,"r")
    ncf_t   = netcdf.netcdf_file(fnm_t,"r")
    ncf_tp1 = netcdf.netcdf_file(fnm_tp1,"r")

    state_tm1 = np.array(ncf_tm1.variables["State"].data) 
    state_t   = np.array(ncf_t.variables["State"].data) 
    state_tp1 = np.array(ncf_tp1.variables["State"].data) 

    obs = open("lorenz96obs_"+str(h).zfill(7)+".txt","w")
    n = 1
    obs.write("%d\n" % len(lobs))
    for l in lobs:
	(x,y)=(l.strip()).split(",")
        if (method > 1):
            if ((n-1)%4==0):
                t=1
                state = state_tm1[int(x)-1]
            elif ((n+1)%4==0):
                t=3
                state = state_tp1[int(x)-1]
            else:
                t=2
                state = state_t[int(x)-1]
        else:
            state = state_t[int(x)-1]
        if (method < 3):
            t=1
        obs.write("%d, %8.3f, %8.3f\n" % (t,float(y),state))
        n=n+1
    obs.close
    ncf_tm1.close
    ncf_t.close
    ncf_tp1.close
 
