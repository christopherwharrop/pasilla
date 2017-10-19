#!/bin/env python
##!//contrib/anaconda/2.3.0/bin/python

vx=open("verify_l96_summary.txt","w")


#for h in range(0,361,36):
#for h in range(0,1201,60):
for h in range(0,1201,60):
  toton=[]
  totof=[]

#  for m in range(1,5):
#  for m in range(1,5):
  for m in (1,2,3):
#  for m in [1]:

    hh=str(h).zfill(4)
    mm=str(m)

    fon=open("vx"+hh+"_"+mm+".txt",'r')
    von=fon.readlines()
    fon.close()
    on=0.0
    for f in von:
	on=on+float(f.split('=')[1]) 
    on=on/float(len(von))
    print m, on
    toton.append(str("%8.3f" % round(float(on),3)))

#    fon=open("vx"+hh+"_"+mm+"_on.txt",'r')
#    von=fon.readlines()
#    fon.close()
#    on=0.0
#    for f in von:
#	on=on+float(f.split('=')[1]) 
#    on=on/float(len(von))
#    toton.append(str("%8.3f" % round(float(on),3)))

#    fof=open("vx"+hh+"_"+mm+"_off.txt",'r')
#    vof=fof.readlines()
#    fof.close()
#    off=0.0
#    for f in vof:
#	off=off+float(f.split('=')[1]) 
#    off=off/float(len(vof))
#    totof.append(str("%8.3f" % round(float(off),3)))

  print str(h)+","+(",".join(toton))
  vx.write(str(h)+","+(",".join(toton))+"\n")
#  print str(h)+","+(",".join(totof))+","+(",".join(toton))
#  vx.write(str(h)+","+(",".join(totof))+","+(",".join(toton))+"\n")
# vx.write("%3s, %8s, %8s, %8s, %8s, %8s, %8s, %8s, %8s \n" % (str(h),(t for t in totof),(t for t in toton)))

vx.close
 
