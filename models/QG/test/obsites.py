#!//contrib/anaconda/2.3.0/bin/python

from scipy.io import netcdf
import numpy as np

text=open("obsites.txt",'r')
lobs=text.readlines()
text.close()

pres=["200", "500", "800"]
assim_window=3

for h in range(36,8639,36):
    fnm_tm1 ="../nature/qgout_"+str(h-assim_window).zfill(7)+".nc"
    fnm_t   ="../nature/qgout_"+str(h).zfill(7)+".nc"
    fnm_tp1 ="../nature/qgout_"+str(h+assim_window).zfill(7)+".nc"

    print h,fnm_tm1, fnm_t, fnm_tp1

    ncf_tm1 = netcdf.netcdf_file(fnm_tm1,"r")
    ncf_t   = netcdf.netcdf_file(fnm_t,"r")
    ncf_tp1 = netcdf.netcdf_file(fnm_tp1,"r")

    psi_tm1 = np.array(ncf_tm1.variables["Psig"].data) 
    psi_t   = np.array(ncf_t.variables["Psig"].data) 
    psi_tp1 = np.array(ncf_tp1.variables["Psig"].data) 

    obs = open("qgobs_"+str(h).zfill(7)+".txt","w")
    n = 1
    obs.write("912 \n")
    for l in lobs:
        if ((n-1)%4==0):
            t=1
        elif ((n+1)%4==0):
            t=3
        else:
            t=2
	(x,y,s,o,a)=(l.strip()).split(",")
        for z in range(2,-1,-1):
             if (t==1):
                 psi = psi_tm1[z,(int(y)-1),(int(x)-1)]
             elif (t==3):
                 psi = psi_tp1[z,(int(y)-1),(int(x)-1)]
             else:
                 psi = psi_t[z,(int(y)-1),(int(x)-1)]
             obs.write("%3s, %8s, %8s, %4s, %14s, %5s \n" % (str("1"),str("%8.3f" % round(float(a),3)),str("%8.3f" % round(float(o),3)),pres[z],str("%13.0f" % psi),s))

        n=n+1
    obs.close
    ncf_tm1.close
    ncf_t.close
    ncf_tp1.close
 
