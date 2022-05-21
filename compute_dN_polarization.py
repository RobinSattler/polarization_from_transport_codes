# 14/05/2022

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

#omega zx array
dozx=0.0125
omin=-0.5
omax=0.4
nozx=int(math.floor((omax-omin)/dozx))

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
dN_dtdozx=np.zeros((nt,nozx),dtype=np.float64)
dN_Ndtdozx=np.zeros((nt,nozx),dtype=np.float64)
hadrons_in_zx=np.zeros(nt,dtype=np.float64)
hadrons_in_tozx=np.zeros(nt,dtype=np.float64)
dN_dxdz=np.zeros((nt,nx,nz),dtype=np.float64)
dN_Ndxdz=np.zeros((nt,nx,nz),dtype=np.float64)
N_hadrons=np.zeros(nt,dtype=np.int64)

if(inputfile[-3:]==".gz"):
  print("Opening gzipped file "+inputfile)
  pi=gzip.open(inputfile,"r")
else:
  print("Opening file "+inputfile)
  pi=open(inputfile,"r")

lines=0
hadrons=0

h_index_for_zx=int(math.floor(15/dt))

for line in pi:
    stuff=line.split() 
    lines=lines+1
    t,x,y,z,pt,rap,sx,sy,sz,ozx=np.float64(stuff[:]) 
    if(t < tmax):
     if((np.abs(rap)<raplim) and (pt>=ptmin) and (pt<=ptmax)):
      hadrons=hadrons+1
      h=int(math.floor(t/dt))
      dN[h]=dN[h]+1
      avgs[h,:]=avgs[h,:]+[sx,sy,sz]
      i=int(math.floor((x+xside)/dx))
      k=int(math.floor((z+zside)/dz))
      m=int(math.floor((ozx-omin)/dozx))
      if((i>=0) and (i<nx) and (k>=0) and (k<nz)):
           dN_dxdz[h,i,k]+=1
           hadrons_in_zx[h]+=1
      if((m>=0) and (m<nozx)):
           dN_dtdozx[h,m]+=1
           hadrons_in_tozx[h]+=1
pi.close()

dN_Ndt=dN/(hadrons*dt)
for h in range(nt):
    if dN[h] != 0:
        avgs[h,:]/=dN[h]
    if hadrons_in_tozx[h] != 0:
        dN_Ndtdozx[h,m] = dN_dtdozx[h,m]/(hadrons_in_tozx[h] * dt * dozx)
    if hadrons_in_zx[h] != 0:
        dN_Ndxdz[h,:,:] = dN_dxdz[h,:,:]/(hadrons_in_zx[h] * dx * dz)

N_hadrons_in_ozx = np.sum(hadrons_in_tozx)
if N_hadrons_in_ozx > 0:
    dN_Ndozx=np.sum(dN_dtdozx,axis=0)/N_hadrons_in_ozx

print("Total Lambdas: "+str(lines)+"\n")
print("Accepted Lambdas: "+str(hadrons)+"\n")

#now we print the results into the output files
sp="          "
cf='{:5.2f}'
of='{:7.4f}'
df='{:8.5e}'

fout=open(outputprefix+"_time_distr.txt","w")
fout.write("#N="+str(hadrons)+"\n")
fout.write("#t      dN_Ndt      dN      <S_x>     <S_y>     <S_z>\n")
for h in range(nt):
    fout.write(cf.format((h+0.5)*dt)+sp+df.format(dN_Ndt[h])+sp+df.format(dN[h])+sp+df.format(avgs[h,0])+sp+df.format(avgs[h,1])+sp+df.format(avgs[h,2])+"\n")
fout.close()


for h in range(nt):
    timestr='{:06.2f}'.format(h*dt+0.5)
    fout=open(outputprefix+"_"+timestr+"_zx_distr.txt","w")
    fout.write("# N = "+str(hadrons_in_zx[h])+"\n")
    fout.write("# time = "+cf.format((h+0.5)*dt)+"\n")
    fout.write("# z      x      dN/(N dx dz)\n")
    for i in range(nx):
        for k in range(nz):
            fout.write(cf.format((i+0.5)*dx)+sp+df.format((k+0.5)*dz)+sp+df.format(dN_Ndxdz[h,i,k])+"\n")
        fout.write("\n")
    fout.close()

fout=open(outputprefix+"_tomegazx_distr.txt","w")
fout.write("# N = "+str(N_hadrons_in_ozx)+"\n")
fout.write("# time = "+cf.format((h+0.5)*dt)+"\n")
fout.write("# t      omega_zx      dN/(N dt domega_zx)\n")
for h in range(nt):
    for m in range(nozx):
        fout.write(cf.format((h+0.5)*dt)+sp+df.format((m+0.5)*dozx)+sp+df.format(dN_Ndtdozx[h,m])+"\n")
    fout.write("\n")
fout.close()


fout=open(outputprefix+"_omegazx_distr.txt","w")
fout.write("# N = "+str(N_hadrons_in_ozx)+"\n")
fout.write("# omega_zx      dN_Ndomegazx\n")
for i in range(nozx):
    fout.write(of.format((i+0.5)*dozx+omin)+sp+df.format(dN_Ndozx[i])+"\n")
fout.close()
