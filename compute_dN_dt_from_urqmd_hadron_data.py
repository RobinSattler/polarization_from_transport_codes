# in this version we assume that we deal with two hadrons with id 27 and 40

import math
import numpy as np
import sys
import os
import gzip


#time resolution
dt=0.5

raplim=1
ptmin=0.1
ptmax=3.

#time max
tmax=160

#number of timesteps (automatically set from dt and tmax)
nt=int(math.floor(tmax/dt))

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files<2):
  print ('Syntax: ./compute_dNpy <inputfile> <output_suffix>')
  sys.exit(1)

inputfile=sys.argv[1]
output_suff=sys.argv[2]

#first, we check that files exist
if(not(os.path.isfile(inputfile))):
  print(inputfile+" does not exist. I quit.\n")
  sys.exit(1)

dN=np.zeros((nt,2),dtype=np.float64)
dN_Ndt=np.zeros((nt,2),dtype=np.float64)

if(inputfile[-3:]==".gz"):
  print("Opening gzipped file "+inputfile)
  pi=gzip.open(inputfile,"r")
else:
  print("Opening file "+inputfile)
  pi=open(inputfile,"r")

lines=0
lambdas=0
sigmas=0
events=0

for line in pi:
    stuff=line.split() 
    lines=lines+1
    #if(not(stuff[0]==b"27") or (stuff[0]==b"40")): (old format)
    if(not(stuff[0]=="27") or (stuff[0]=="40")): 
        if stuff[0] == "event":
            events+=1
        continue
    t,x,y,z,En,px,py,pz=np.float64(stuff[1:]) #here px, py and pz is the average spin, not the momentum
    if(t < tmax):
      rapidity=0.5*math.log((En+pz)/(En-pz))
      pt=math.sqrt(px**2+py**2)
      if((abs(rapidity)<raplim) and (pt>=ptmin) and (pt<=ptmax)):
         h=int(math.floor(t/dt))
         #if(stuff[0]==b"27"): (old format)
         if(stuff[0]=="27"):
            lambdas=lambdas+1
            dN[h,0]=dN[h,0]+1
         else:   
            sigmas=sigmas+1
            dN[h,1]=dN[h,0]+1

pi.close()

if(lambdas > 0):
   dN_Ndt[:,0]=dN[:,0]/(lambdas*dt)
if(sigmas > 0):
   dN_Ndt[:,1]=dN[:,1]/(sigmas*dt)

#now we print the results into a file
sp="          "
fout=open("Lambda_time_distr_"+output_suff+".dat","w")
fout.write("#N Lambdas = "+str(lambdas)+"\n")
fout.write("#N events = "+str(events)+"\n")
fout.write("#t      dN/(N_Lambdas dt)      dN/(N_events dt)\n")
for h in range(nt):
    fout.write('{:5.2f}'.format((h+0.5)*dt)+sp+'{:7.4e}'.format(dN_Ndt[h,0])+sp+'{:7.4e}'.format(dN[h,0]/(dt*events))+"\n")
fout.close()

fout=open("Sigma_time_distr_"+output_suff+".dat","w")
fout.write("#N Sigmas = "+str(sigmas)+"\n")
fout.write("#N events = "+str(events)+"\n")
fout.write("#t      dN/(N_Sigmas dt)      dN/(N_events dt)\n")
for h in range(nt):
    fout.write('{:5.2f}'.format((h+0.5)*dt)+sp+'{:7.4e}'.format(dN_Ndt[h,1])+sp+'{:7.4e}'.format(dN[h,1]/(dt*events))+"\n")
fout.close()
