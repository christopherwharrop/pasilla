#!/bin/env python

import sys
from scipy.io import netcdf
import numpy as np

#method=int(sys.argv[1])
#print "Generating obs for method = ", method

text=open("obs_sites.txt",'r')
lobs=text.readlines()
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

    # Pre-calculate random errors from N(0,1) to apply to the observations
    # ensuring that each method has the same obs values
    err = np.random.normal(0,1,len(lobs))

    for method in (1,2,3):

        obs = open("lorenz96obs_"+str(method)+"_"+str(h).zfill(7)+".txt","w")
        if (method > 1):
            obs.write("%d\n" % (3 * len(lobs)))
        else:
            obs.write("%d\n" % len(lobs))
        for l in lobs:
            (x,y)=(l.strip()).split(",")
            if (method > 1):
                for t in (1,2,3):
                    if (t==1):
                        state = state_tm1[int(x)-1]
                    elif (t==3):
                        state = state_tp1[int(x)-1]
                    else:
                        state = state_t[int(x)-1]
                    if (method==2):
#                        obs.write("%d, %8.3f, %8.3f\n" % (1,float(y),state))
                        obs.write("%d, %f, %20.15f\n" % (1,float(y),state))
                    else:
#                        obs.write("%d, %8.3f, %8.3f\n" % (t,float(y),state))
                        obs.write("%d, %f, %20.15f\n" % (t,float(y),state))

            else:
                state = state_t[int(x)-1]
#                obs.write("%d, %8.3f, %8.3f\n" % (1,float(y),state))
                obs.write("%d, %f, %20.15f\n" % (1,float(y),state))
        obs.close

    ncf_tm1.close
    ncf_t.close
    ncf_tp1.close
 
