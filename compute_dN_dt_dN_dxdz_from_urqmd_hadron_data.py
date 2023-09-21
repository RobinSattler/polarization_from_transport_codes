#!/usr/bin/env python3

import math
import numpy as np
import sys
import os
import gzip

# hadron ids
hadron_ids=(27,40)
hadron_names=("Lambda","Sigma")

# time resolution
dt=1
# time start
time_min=0.5

# x,z resolution
dx=1
dz=1

# histograms go from -xside to xside
xside=20
zside=20

raplim=100
ptmin=0.1
ptmax=3.

# time max at center cell
tmax=160

# preliminary consistency check
if len(hadron_ids) != len(hadron_names):
    print("Error, mismatch between hadron_ids and hadron_names tuples. Please, check the script!")
    sys.exit(1)
else:
    nhad=len(hadron_ids)

#number of timesteps (automatically set from dt and tmax)
nt=int(math.floor((tmax+dt/2.-time_min)/dt))

nx=int(math.floor(2*xside/dx))
nz=int(math.floor(2*zside/dz))

# we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files<2):
  print ('Syntax: python3 compute_dN_dt_dN_dxdz_from_urqmd_hadron_data.py <inputfile> <output_suffix>')
  sys.exit(1)

inputfile=sys.argv[1]
output_suff=sys.argv[2]

# first, we check that files exist
if(not(os.path.isfile(inputfile))):
  print(inputfile+" does not exist. I quit.\n")
  sys.exit(1)

dN=np.zeros((nhad,nt),dtype=np.float64)
dN_Ndt=np.zeros((nhad,nt),dtype=np.float64)
hadrons_in_xz=np.zeros((nhad,nt),dtype=np.float64)
dN_dxdz=np.zeros((nhad,nt,nx,nz),dtype=np.float64)
dN_Ndxdz=np.zeros((nhad,nt,nx,nz),dtype=np.float64)
N_hadrons=np.zeros(nhad,dtype=np.int64)

if(inputfile[-3:]==".gz"):
  print("Opening gzipped file "+inputfile)
  pi=gzip.open(inputfile,"rb")
  gzip_data=True
else:
  print("Opening file "+inputfile)
  pi=open(inputfile,"r")
  gzip_data=False

lines=0
events=0

for line in pi:
    stuff=line.split() 
    lines=lines+1
    if gzip_data:
        entry = stuff[0].decode('utf-8')
    else:
        entry = stuff[0]
    if entry == "event":
        events+=1
        continue
    else:
        pid = int(entry)
    #if(not(stuff[0]==b"27") or (stuff[0]==b"40")): (old format)
    if(not(pid in hadron_ids)): 
        continue
    t,x,y,z,En,px,py,pz=np.float64(stuff[1:]) #here px, py and pz is the average spin, not the momentum
    if(t < tmax):
      rapidity=0.5*math.log((En+pz)/(En-pz))
      pt=math.sqrt(px**2+py**2)
      if((abs(rapidity)<raplim) and (pt>=ptmin) and (pt<=ptmax)):
         h=int(math.floor((t-time_min)/dt))
         #if(stuff[0]==b"27"): (old format)
         i=int(math.floor((x+xside)/dx))
         k=int(math.floor((z+zside)/dz))
         indx_had=hadron_ids.index(pid)
         N_hadrons[indx_had]+=1
         dN[indx_had,h]=dN[indx_had,h]+1
         if((i>=0) and (i<nx) and (k>=0) and (k<nz)):
             dN_dxdz[indx_had,h,i,k]+=1
             hadrons_in_xz[indx_had,h]+=1

pi.close()

for i in range(nhad):
    if(N_hadrons[i] > 0):
        dN_Ndt[i,:]=dN[i,:]/(N_hadrons[i]*dt)
    for h in range(nt):
        if hadrons_in_xz[i,h] > 0:
            dN_Ndxdz[i,h,:,:]=dN_dxdz[i,h,:,:]/(hadrons_in_xz[i,h]*dx*dz*dt)

#now we print the results into a file
sp="          "
cf='{:5.2f}'
df='{:8.5e}'
for q in range(nhad):
    hadname=hadron_names[q]
    fout=open(hadname+"_time_distr_"+output_suff+".dat","w")
    fout.write("# N "+hadname+" = "+str(N_hadrons[q])+"\n")
    fout.write("# N events = "+str(events)+"\n")
    fout.write("# t      dN/(N_"+hadname+" dt)      dN/(N_events dt)\n")
    for h in range(nt):
        fout.write(cf.format(h*dt+0.5+time_min)+sp+df.format(dN_Ndt[0,h])+sp+df.format(dN[0,h]/(dt*events))+"\n")
    fout.close()
    for h in range(nt):
        timestr='{:06.2f}'.format(h*dt+0.5)
        fout=open(hadname+"_zx_distr_t_"+timestr+"_"+output_suff+".dat","w")
        fout.write("# N "+hadname+" = "+str(hadrons_in_xz[q,h])+"\n")
        fout.write("# N events = "+str(events)+"\n")
        fout.write("# time = "+cf.format((h+0.5)*dt)+"\n")
        fout.write("# z    x      dN/(N_"+hadname+" dxdzdt)      dN/(N_events dxdzdt)\n")
        for i in range(nx):
            for k in range(nz):
                fout.write(cf.format((k+0.5)*dz-zside)+sp+cf.format((i+0.5)*dx-xside)+sp+df.format(dN_Ndxdz[q,h,i,k])+sp+df.format(dN_dxdz[q,h,i,k]/(events*dx*dz*dt))+"\n")
            fout.write("\n")
        fout.close()
