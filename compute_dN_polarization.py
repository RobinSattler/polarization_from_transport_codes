# 10/05/2022

import math
import numpy as np
import sys
import os
import gzip

#time resolution
dt=0.5

#time max
tmax=100

#rapidity limit in absolute value
raplim=1.0

#pt range
ptmin=0.1
ptmax=3

#x,z resolution
dx=0.5
dz=0.5

#histograms go from -xside to xside
xside=20
zside=20

#number of timesteps (automatically set from dt and tmax)
nt=int(math.floor(tmax/dt))

nx=int(math.floor(2*xside/dx))
nz=int(math.floor(2*zside/dz))

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files<2):
  print ('Syntax: ./compute_dNpy <inputfile> <output_prefix>')
  sys.exit(1)

inputfile=sys.argv[1]
outputprefix=sys.argv[2]

#first, we check that files exist
if(not(os.path.isfile(inputfile))):
  print(inputfile+" does not exist. I quit.\n")
  sys.exit(1)

dN=np.zeros(nt,dtype=np.float64)
avgs=np.zeros((nt,3),dtype=np.float64)
dN_dxdz=np.zeros((nx,nz),dtype=np.float64)

if(inputfile[-3:]==".gz"):
  print("Opening gzipped file "+inputfile)
  pi=gzip.open(inputfile,"r")
else:
  print("Opening file "+inputfile)
  pi=open(inputfile,"r")

lines=0
hadrons=0
hadrons_in_xz=0

h_index_for_xz=int(math.floor(15/dt))

for line in pi:
    stuff=line.split() 
    lines=lines+1
    t,x,y,z,pt,rap,sx,sy,sz=np.float64(stuff[:]) 
    if(t < tmax):
     if((np.abs(rap)<raplim) and (pt>=ptmin) and (pt<=ptmax)):
      hadrons=hadrons+1
      h=int(math.floor(t/dt))
      dN[h]=dN[h]+1
      avgs[h,:]=avgs[h,:]+[sx,sy,sz]
      i=int(math.floor((x+xside)/dx))
      k=int(math.floor((z+zside)/dz))
      if((h==h_index_for_xz) and (i>=0) and (i<nx) and (k>=0) and (k<nz)):
             dN_dxdz[i,k]=dN_dxdz[i,k]+1
             hadrons_in_xz=hadrons_in_xz+1

pi.close()

dN_Ndt=dN/(hadrons*dt)
dN_Ndxdz=dN_dxdz/(hadrons_in_xz*dx*dz)
for i in range(nt):
    if dN[i] != 0:
        avgs[i,:]=avgs[i,:]/dN[i]

print("Total Lambdas: "+str(lines)+"\n")
print("Accepted Lambdas: "+str(hadrons)+"\n")

#now we print the results into a file
sp="          "
fout=open(outputprefix+"_time_distr.txt","w")
fout.write("#N="+str(hadrons)+"\n")
fout.write("#t      dN_Ndt      dN      <S_x>     <S_y>     <S_z>\n")
for h in range(nt):
    fout.write('{:5.2f}'.format((h+0.5)*dt)+sp+'{:7.4f}'.format(dN_Ndt[h])+sp+'{:7.4f}'.format(dN[h])+sp+'{:7.5e}'.format(avgs[h,0])+sp+'{:7.5e}'.format(avgs[h,1])+sp+'{:7.5e}'.format(avgs[h,2])+"\n")
fout.close()

fout=open(outputprefix+"_xz_distr.txt","w")
fout.write("#N="+str(hadrons_in_xz)+"\n")
fout.write("#x      z      dN_Ndxdz\n")
for k in range(nz):
    for i in range(nx):
        fout.write('{:5.2f}'.format((i+0.5)*dx)+sp+'{:5.2f}'.format((k+0.5)*dz)+sp+'{:7.4f}'.format(dN_Ndxdz[i,k])+"\n")
    fout.write("\n")
fout.close()
