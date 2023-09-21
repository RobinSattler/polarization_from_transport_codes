# version 0.2.0 - 22/04/2022

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import gzip
from scipy import interpolate


"""
   It computes the vorticity from the data produced by SMASH with thermodynamic lattice output,
   rho_eckart, tmn_landau, landau_velocity, j_QBS, also if gzipped.
   A tabulated EoS in ascii format with mapping
   energy density, baryon density and electric charge density to
   temperature, pressure, baryon, charge and strangeness chemical potentials
   must be provided, as well.
   This version does not use numpy to compute partial derivatives.
"""

# PARAMETERS
eos_type = 2 #1 = Official SMASH HRG (e,nB,nQ), 2 = UrMQD

eos_dir = "EOS_HG_UrQMD"

der_type = 2 #1 = numpy derivatives, 2 = 2nd order centered differences with 1st at borders

verbose = True # if True it informs about the advancement of the program

temp_limit = 0.001 # temperature limit to accept a cell in GeV

#we set the parameter hbarc
hbarc=0.197326

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if (N_input_files!=2):
   print ('Syntax: ./compute_vorticity_myder_smash.py <pickled archive file with SMASH data> <outputfile>')
   sys.exit(1)

inputfile=sys.argv[1]
outputfile=sys.argv[2]

if (inputfile[-3:]==".gz"):
    if (verbose):
       print("Opening gzipped file "+inputfile)
    pi=gzip.open(inputfile,"rb")
else:
    if (verbose):
       print("Opening file "+inputfile)
    pi=open(inputfile,"rb")

indata=pickle.load(pi)
pi.close()
lattice,tt,N_events,Tmunu,J,v = indata[:]


dt=tt[1]-tt[0]
nt=len(tt)
nx,ny,nz=lattice["dimensions"]
dx,dy,dz=lattice["spacing"]
xstart,ystart,zstart=lattice["origin"]

xx=np.linspace(xstart+dx/2,xstart+dx*(nx-0.5),nx)
yy=np.linspace(ystart+dy/2,ystart+dy*(ny-0.5),ny)
zz=np.linspace(zstart+dz/2,zstart+dz*(nz-0.5),nz)

if (verbose):
    print("Coarse graining data read, now diving by the number of events, i.e.: "+str(N_events))
Tmunu=Tmunu/N_events
J=J/N_events
v=v/N_events

# we use these arrays for convenience
vx=v[:,0,:,:,:]
vy=v[:,1,:,:,:]
vz=v[:,2,:,:,:]

# important indexes

# energy density in T_munu
iT00=0

#indexes corresponding to jQ[0:3] in j_QBS
kQ0=0
kQ1=1
kQ2=2
kQ3=3
#indexes corresponding to jB[0:3] in j_QBS
kB0=4
kB1=5
kB2=6
kB3=7
#indexes corresponding to jS[0:3] in j_QBS
kS0=8
kS1=9
kS2=10
kS3=11

if (verbose):
    print("Computing the gamma Lorentz factor")
glf=1/np.sqrt(1-np.sum(v**2,axis=1))

if (verbose):
    print("Computing rhoQ")
rhoQ=glf*(J[:,kQ0,:,:,:]-J[:,kQ1,:,:,:]*vx-J[:,kQ2,:,:,:]*vy-J[:,kQ3,:,:,:]*vz)

if (verbose):
    print("Computing rhoB")
rhoB=glf*(J[:,kB0,:,:,:]-J[:,kB1,:,:,:]*vx-J[:,kB2,:,:,:]*vy-J[:,kB3,:,:,:]*vz)

if (verbose):
    print("Computing rhoS")
rhoS=glf*(J[:,kS0,:,:,:]-J[:,kS1,:,:,:]*vx-J[:,kS2,:,:,:]*vy-J[:,kS3,:,:,:]*vz)



if eos_type == 1: #smash eos type

    eosfile=eos_dir+"/hadgas_eos_SMASH.dat"

    en_list=[]
    B_list=[]
    Q_list=[]

    def count_points(eosfile):
        if(eosfile[-3:]==".gz"):
            if(verbose):
                print("Opening gzipped EoS file "+eosfile)
            fin=gzip.open(eosfile,"rb")
        else:
            if(verbose):
                print("Opening EoS file "+eosfile)
            fin=open(eosfile,"rb")
        # we skip the first row
        fin.readline()
        # we read the first row
        eref,bref,qref=np.float64(fin.readline().split()[0:3])
        en_list.append(eref)
        B_list.append(bref)
        Q_list.append(qref)
        nE=1;nB=1;nQ=1
        for line in fin:
            enval,bval,qval=np.float64(line.split()[0:3])
            if ((nE==1) and (nB==1)):
                if (bval==bref):
                    nQ=nQ+1
                    Q_list.append(qval)
                else:
                    nB=nB+1
                    B_list.append(bval)
                    bref=bval
            if ((nE==1) and (enval==eref)):
                if (bval!=bref):
                    nB=nB+1
                    B_list.append(bval)
                    bref=bval
            if (enval!=eref):
                nE=nE+1
                en_list.append(enval)
                eref=enval
        fin.close()
        return nE,nB,nQ

    if (verbose):
        print("Retrieving information about the tabulated EoS structure")
    nE,nB,nQ=count_points(eosfile)

    en_arr=np.array(en_list,dtype=np.float64)
    B_arr=np.array(B_list,dtype=np.float64)
    Q_arr=np.array(Q_list,dtype=np.float64)

    if (verbose):
        print("Retrieving the data from the tabulated EoS")
    raw_data=np.loadtxt(eosfile,skiprows=1,usecols=(3,4,5,6,7))
    data=raw_data.reshape((nE,nB,nQ,5))

    if (verbose):
        print("Computing the temperature by interpolating the EoS")
    temp = interpolate.interpn((en_arr,B_arr,Q_arr),data[:,:,:,0],(Tmunu[:,iT00,:,:,:],rhoB,rhoQ),bounds_error=False,fill_value=-1)

elif eos_type == 2: #urqmd eos type 

    #these are the unit values of energy density and pressure
    e0=0.14651751415742
    #these are the unit values of the net baryon density and entropy density
    n0=0.15891

    fstd=eos_dir+"/hadgas_eos.dat"
    Ne_std=2001
    Nn_std=401
    en_std_max=1000.
    rho_std_max=40.
    enarr_std=np.linspace(0.,en_std_max,num=Ne_std)
    rhoarr_std=np.linspace(0.,rho_std_max,num=Nn_std)

    fmed=eos_dir+"/hg_eos_small.dat"
    Ne_med=201
    Nn_med=201
    en_med_max=10.
    rho_med_max=2.
    enarr_med=np.linspace(0.,en_med_max,num=Ne_med)
    rhoarr_med=np.linspace(0.,rho_med_max,num=Nn_med)

    fmin=eos_dir+"/hg_eos_mini.dat"
    Ne_min=201
    Nn_min=201
    en_min_max=0.1
    rho_min_max=0.02
    enarr_min=np.linspace(0.,en_min_max,num=Ne_min)
    rhoarr_min=np.linspace(0.,rho_min_max,num=Nn_min)

    temparr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)
    muarr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)
    sarr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)
    parr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)

    temparr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)
    muarr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)
    sarr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)
    parr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)

    temparr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)
    muarr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)
    sarr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)
    parr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)

    def readeos(ff,tarr,marr,parr,sarr,nx,ny):
        for j in range(ny):
            for i in range(nx):
                stuff=ff.readline().split()
                tarr[i,j],marr[i,j],parr[i,j],sarr[i,j]=np.float64(stuff[0]),np.float64(stuff[1]),np.float64(stuff[3]),np.float64(stuff[5])

    if verbose:
        print("Reading the tabulated EoS from the files")

    with open(fstd,"r") as infile:
        readeos(infile,temparr_std,muarr_std,parr_std,sarr_std,Ne_std,Nn_std)
        temp_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, temparr_std.transpose(), kind='linear')
        muB_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, muarr_std.transpose(), kind='linear')
        p_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, parr_std.transpose(), kind='linear')
        s_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, sarr_std.transpose(), kind='linear')

    with open(fmed,"r") as infile:
        readeos(infile,temparr_med,muarr_med,parr_med,sarr_med,Ne_med,Nn_med)
        temp_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, temparr_med.transpose(), kind='linear')
        muB_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, muarr_med.transpose(), kind='linear')
        p_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, parr_med.transpose(), kind='linear')
        s_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, sarr_med.transpose(), kind='linear')

    with open(fmin,"r") as infile:
        readeos(infile,temparr_min,muarr_min,parr_min,sarr_min,Ne_min,Nn_min)
        temp_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, temparr_min.transpose(), kind='linear')
        muB_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, muarr_min.transpose(), kind='linear')
        p_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, parr_min.transpose(), kind='linear')
        s_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, sarr_min.transpose(), kind='linear')
 
    if verbose:
        print("Done.\n")

    def get_mub_T(rhoB_input_w_sign,edens):
        #before callin this function we already checked that both arguments are > 0
        compute=True
        if rhoB_input_w_sign >= 0:
            rhoB_tmp=rhoB_input_w_sign
        else:
            rhoB_tmp = 0.
        edens_val=edens/e0
        rhoB_val=rhoB_tmp/n0
        if(edens_val<=en_std_max):
            if((edens_val<en_min_max ) and (rhoB_val<rho_min_max)):
                ftemp=temp_interp_min
                fmuB=muB_interp_min
                fpr=p_interp_min
                fs=s_interp_min
            if(edens_val<en_med_max ) and (rhoB_val<rho_med_max) and ((edens_val>=en_min_max ) or (rhoB_val>=rho_min_max)):     
                ftemp=temp_interp_med
                fmuB=muB_interp_med
                fpr=p_interp_med
                fs=s_interp_med
            if((edens_val>=en_med_max ) or (rhoB_val>=rho_med_max)):     
                if(rhoB_val>rho_std_max):
                    print("Net baryon density exceeding the maximum of the table. Changed from "+str(rhoB_val)+" to "+str(rho_std_max*0.999999))
                    rhoB_val=rho_std_max*0.999999
                ftemp=temp_interp_std
                fmuB=muB_interp_std
                fpr=p_interp_std
                fs=s_interp_std
        elif(edens_val>en_std_max):    
            compute=False
#           temperature=350./1000.
#           muB=3./1000.
            temperature=0.
            muB=0.
            pressure=0.
            entr_dens=0.
        else:        
            compute=False
            temperature=0.
            muB=0.
            pressure=0.
            entr_dens=0.
     
        if(compute): #all is expressed in GeV
            temperature=ftemp(edens_val,rhoB_val)[0]/1000.
            #muB=3*fmuB(edens,rhoB)[0]/1000. 
            #pressure=fpr(edens,rhoB)[0]*e0
            #entr_dens=fs(edens,rhoB)[0]*n0

        #return muB, temperature, pressure, entr_dens
        return temperature

    temp=np.zeros((nt,nx,ny,nz),dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    temp[h,i,j,k]=get_mub_T(rhoB[h,i,j,k],Tmunu[h,iT00,i,j,k])


bt=np.zeros(temp.shape,dtype=np.float64)
bx=np.zeros(temp.shape,dtype=np.float64)
by=np.zeros(temp.shape,dtype=np.float64)
bz=np.zeros(temp.shape,dtype=np.float64)

fout=open("cells.dat","w")
sp="     "
for h in range(nt):
    count_null=0
    count_neg=0
    count_pos=0
    count_good=0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if temp[h,i,j,k] < 0:
                    count_neg=count_neg+1
                    continue
                if temp[h,i,j,k] == 0:
                    count_null=count_null+1
                    continue
                if temp[h,i,j,k] > 0:
                    #print(str(Tmunu[h,iT00,i,j,k])+sp+str(rhoB[h,i,j,k])+sp+str(rhoQ[h,i,j,k])+sp+str(temp[h,i,j,k]))
                    count_pos=count_pos+1
                if(temp[h,i,j,k]>temp_limit):
                    count_good=count_good+1
                    bt[h,i,j,k]=glf[h,i,j,k]/temp[h,i,j,k]
                    #we are using the covariant components, index down, so they get a -1 sign (Minkowski metric signature +---)
                    bx[h,i,j,k]=-vx[h,i,j,k]*glf[h,i,j,k]/temp[h,i,j,k]
                    by[h,i,j,k]=-vy[h,i,j,k]*glf[h,i,j,k]/temp[h,i,j,k]
                    bz[h,i,j,k]=-vz[h,i,j,k]*glf[h,i,j,k]/temp[h,i,j,k]


    fout.write(str(tt[h])+sp+str(count_neg)+sp+str(count_null)+sp+str(count_pos)+sp+str(count_good)+"\n")
    
fout.close()


if der_type == 1:

    dbt_dx=np.gradient(bt,dx,axis=1)
    dbt_dy=np.gradient(bt,dy,axis=2)
    dbt_dz=np.gradient(bt,dz,axis=3)
    dbx_dt=np.gradient(bx,dt,axis=0)
    dby_dt=np.gradient(by,dt,axis=0)
    dbz_dt=np.gradient(bz,dt,axis=0)

    dbx_dy=np.gradient(bx,dy,axis=2)
    dbx_dz=np.gradient(bx,dz,axis=3)
    dby_dx=np.gradient(by,dx,axis=1)
    dby_dz=np.gradient(by,dz,axis=3)
    dbz_dx=np.gradient(bz,dx,axis=1)
    dbz_dy=np.gradient(bz,dy,axis=2)

elif der_type == 2:

    dbt_dx=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(i==0):
                        if((bt[h,0,j,k]!=0) and (bt[h,1,j,k]!=0)):
                            dbt_dx[h,0,j,k]=(bt[h,1,j,k]-bt[h,0,j,k])/dx
                    elif(i==nx-1):
                        if((bt[h,nx-1,j,k]!=0) and (bt[h,nx-2,j,k]!=0)):
                            dbt_dx[h,nx-1,j,k]=(bt[h,nx-1,j,k]-bt[h,nx-2,j,k])/dx
                    elif((bt[h,i-1,j,k]!=0) and (bt[h,i+1,j,k]!=0)):
                        dbt_dx[h,i,j,k]=(bt[h,i+1,j,k]-bt[h,i-1,j,k])/(2*dx)
                    elif((bt[h,i-1,j,k]!=0) and (bt[h,i,j,k]!=0)):
                        dbt_dx[h,i,j,k]=(bt[h,i,j,k]-bt[h,i-1,j,k])/(dx)
                    elif((bt[h,i+1,j,k]!=0) and (bt[h,i,j,k]!=0)):
                        dbt_dx[h,i,j,k]=(bt[h,i+1,j,k]-bt[h,i,j,k])/(dx)

    dbt_dy=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(j==0):
                        if((bt[h,i,0,k]!=0) and (bt[h,i,1,k]!=0)):
                            dbt_dy[h,i,0,k]=(bt[h,i,1,k]-bt[h,i,0,k])/dy
                    elif(j==ny-1):
                        if((bt[h,i,ny-1,k]!=0) and (bt[h,i,ny-2,k]!=0)):
                            dbt_dy[h,i,ny-1,k]=(bt[h,i,ny-1,k]-bt[h,i,ny-2,k])/dy
                    elif((bt[h,i,j-1,k]!=0) and (bt[h,i,j+1,k]!=0)):
                        dbt_dy[h,i,j,k]=(bt[h,i,j+1,k]-bt[h,i,j-1,k])/(2*dy)
                    elif((bt[h,i,j-1,k]!=0) and (bt[h,i,j,k]!=0)):
                        dbt_dy[h,i,j,k]=(bt[h,i,j,k]-bt[h,i,j-1,k])/(dy)
                    elif((bt[h,i,j+1,k]!=0) and (bt[h,i,j,k]!=0)):
                        dbt_dy[h,i,j,k]=(bt[h,i,j+1,k]-bt[h,i,j,k])/(dy)

    dbt_dz=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(k==0):
                        if((bt[h,i,j,0]!=0) and (bt[h,i,j,1]!=0)):
                            dbt_dz[h,i,j,0]=(bt[h,i,j,1]-bt[h,i,j,0])/dz
                    elif(k==nz-1):
                        if((bt[h,i,j,nz-1]!=0) and (bt[h,i,j,nz-2]!=0)):
                            dbt_dz[h,i,j,nz-1]=(bt[h,i,j,nz-1]-bt[h,i,j,nz-2])/dz
                    elif((bt[h,i,j,k-1]!=0) and (bt[h,i,j,k+1]!=0)):
                        dbt_dz[h,i,j,k]=(bt[h,i,j,k+1]-bt[h,i,j,k-1])/(2*dz)
                    elif((bt[h,i,j,k-1]!=0) and (bt[h,i,j,k]!=0)):
                        dbt_dz[h,i,j,k]=(bt[h,i,j,k]-bt[h,i,j,k-1])/(dz)
                    elif((bt[h,i,j,k+1]!=0) and (bt[h,i,j,k]!=0)):
                        dbt_dz[h,i,j,k]=(bt[h,i,j,k+1]-bt[h,i,j,k])/(dz)

    dbx_dt=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(h==0):
                        if((bx[0,i,j,k]!=0) and (bx[1,i,j,k]!=0)):
                            dbx_dt[0,i,j,k]=(bx[1,i,j,k]-bx[0,i,j,k])/dt
                    elif(h==nt-1):
                        if((bx[nt-1,i,j,k]!=0) and (bx[nt-2,i,j,k]!=0)):
                            dbx_dt[nt-1,i,j,k]=(bx[nt-1,i,j,k]-bx[nt-2,i,j,k])/dt
                    elif((bx[h-1,i,j,k]!=0) and (bx[h+1,i,j,k]!=0)):
                        dbx_dt[h,i,j,k]=(bx[h+1,i,j,k]-bx[h-1,i,j,k])/(2*dt)
                    elif((bx[h-1,i,j,k]!=0) and (bx[h,i,j,k]!=0)):
                        dbx_dt[h,i,j,k]=(bx[h,i,j,k]-bx[h-1,i,j,k])/(dt)
                    elif((bx[h+1,i,j,k]!=0) and (bx[h,i,j,k]!=0)):
                        dbx_dt[h,i,j,k]=(bx[h+1,i,j,k]-bx[h,i,j,k])/(dt)


    dby_dt=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(h==0):
                        if((by[0,i,j,k]!=0) and (by[1,i,j,k]!=0)):
                            dby_dt[0,i,j,k]=(by[1,i,j,k]-by[0,i,j,k])/dt
                    elif(h==nt-1):
                        if((by[nt-1,i,j,k]!=0) and (by[nt-2,i,j,k]!=0)):
                            dby_dt[nt-1,i,j,k]=(by[nt-1,i,j,k]-by[nt-2,i,j,k])/dt
                    elif((by[h-1,i,j,k]!=0) and (by[h+1,i,j,k]!=0)):
                        dby_dt[h,i,j,k]=(by[h+1,i,j,k]-by[h-1,i,j,k])/(2*dt)
                    elif((by[h-1,i,j,k]!=0) and (by[h,i,j,k]!=0)):
                        dby_dt[h,i,j,k]=(by[h,i,j,k]-by[h-1,i,j,k])/(dt)
                    elif((by[h+1,i,j,k]!=0) and (by[h,i,j,k]!=0)):
                        dby_dt[h,i,j,k]=(by[h+1,i,j,k]-by[h,i,j,k])/(dt)


    dbz_dt=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(h==0):
                        if((bz[0,i,j,k]!=0) and (bz[1,i,j,k]!=0)):
                            dbz_dt[0,i,j,k]=(bz[1,i,j,k]-bz[0,i,j,k])/dt
                    elif(h==nt-1):
                        if((bz[nt-1,i,j,k]!=0) and (bz[nt-2,i,j,k]!=0)):
                            dbz_dt[nt-1,i,j,k]=(bz[nt-1,i,j,k]-bz[nt-2,i,j,k])/dt
                    elif((bz[h-1,i,j,k]!=0) and (bz[h+1,i,j,k]!=0)):
                        dbz_dt[h,i,j,k]=(bz[h+1,i,j,k]-bz[h-1,i,j,k])/(2*dt)
                    elif((bz[h-1,i,j,k]!=0) and (bz[h,i,j,k]!=0)):
                        dbz_dt[h,i,j,k]=(bz[h,i,j,k]-bz[h-1,i,j,k])/(dt)
                    elif((bz[h+1,i,j,k]!=0) and (bz[h,i,j,k]!=0)):
                        dbz_dt[h,i,j,k]=(bz[h+1,i,j,k]-bz[h,i,j,k])/(dt)


    dbx_dy=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(j==0):
                        if((bx[h,i,0,k]!=0) and (bx[h,i,1,k]!=0)):
                            dbx_dy[h,i,0,k]=(bx[h,i,1,k]-bx[h,i,0,k])/dy
                    elif(j==ny-1):
                        if((bx[h,i,ny-1,k]!=0) and (bx[h,i,ny-2,k]!=0)):
                            dbx_dy[h,i,ny-1,k]=(bx[h,i,ny-1,k]-bx[h,i,ny-2,k])/dy
                    elif((bx[h,i,j-1,k]!=0) and (bx[h,i,j+1,k]!=0)):
                        dbx_dy[h,i,j,k]=(bx[h,i,j+1,k]-bx[h,i,j-1,k])/(2*dy)
                    elif((bx[h,i,j-1,k]!=0) and (bx[h,i,j,k]!=0)):
                        dbx_dy[h,i,j,k]=(bx[h,i,j,k]-bx[h,i,j-1,k])/(dy)
                    elif((bx[h,i,j+1,k]!=0) and (bx[h,i,j,k]!=0)):
                        dbx_dy[h,i,j,k]=(bx[h,i,j+1,k]-bx[h,i,j,k])/(dy)


    dbz_dy=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(j==0):
                        if((bz[h,i,0,k]!=0) and (bz[h,i,1,k]!=0)):
                            dbz_dy[h,i,0,k]=(bz[h,i,1,k]-bz[h,i,0,k])/dy
                    elif(j==ny-1):
                        if((bz[h,i,ny-1,k]!=0) and (bz[h,i,ny-2,k]!=0)):
                            dbz_dy[h,i,ny-1,k]=(bz[h,i,ny-1,k]-bz[h,i,ny-2,k])/dy
                    elif((bz[h,i,j-1,k]!=0) and (bz[h,i,j+1,k]!=0)):
                        dbz_dy[h,i,j,k]=(bz[h,i,j+1,k]-bz[h,i,j-1,k])/(2*dy)
                    elif((bz[h,i,j-1,k]!=0) and (bz[h,i,j,k]!=0)):
                        dbz_dy[h,i,j,k]=(bz[h,i,j,k]-bz[h,i,j-1,k])/(dy)
                    elif((bz[h,i,j+1,k]!=0) and (bz[h,i,j,k]!=0)):
                        dbz_dy[h,i,j,k]=(bz[h,i,j+1,k]-bz[h,i,j,k])/(dy)


    dby_dx=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(i==0):
                        if((by[h,0,j,k]!=0) and (by[h,1,j,k]!=0)):
                            dby_dx[h,0,j,k]=(by[h,1,j,k]-by[h,0,j,k])/dx
                    elif(i==nx-1):
                        if((by[h,nx-1,j,k]!=0) and (by[h,nx-2,j,k]!=0)):
                            dby_dx[h,nx-1,j,k]=(by[h,nx-1,j,k]-by[h,nx-2,j,k])/dx
                    elif((by[h,i-1,j,k]!=0) and (by[h,i+1,j,k]!=0)):
                        dby_dx[h,i,j,k]=(by[h,i+1,j,k]-by[h,i-1,j,k])/(2*dx)
                    elif((by[h,i-1,j,k]!=0) and (by[h,i,j,k]!=0)):
                        dby_dx[h,i,j,k]=(by[h,i,j,k]-by[h,i-1,j,k])/(dx)
                    elif((by[h,i+1,j,k]!=0) and (by[h,i,j,k]!=0)):
                        dby_dx[h,i,j,k]=(by[h,i+1,j,k]-by[h,i,j,k])/(dx)

    dbz_dx=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(i==0):
                        if((bz[h,0,j,k]!=0) and (bz[h,1,j,k]!=0)):
                            dbz_dx[h,0,j,k]=(bz[h,1,j,k]-bz[h,0,j,k])/dx
                    elif(i==nx-1):
                        if((bz[h,nx-1,j,k]!=0) and (bz[h,nx-2,j,k]!=0)):
                            dbz_dx[h,nx-1,j,k]=(bz[h,nx-1,j,k]-bz[h,nx-2,j,k])/dx
                    elif((bz[h,i-1,j,k]!=0) and (bz[h,i+1,j,k]!=0)):
                        dbz_dx[h,i,j,k]=(bz[h,i+1,j,k]-bz[h,i-1,j,k])/(2*dx)
                    elif((bz[h,i-1,j,k]!=0) and (bz[h,i,j,k]!=0)):
                        dbz_dx[h,i,j,k]=(bz[h,i,j,k]-bz[h,i-1,j,k])/(dx)
                    elif((bz[h,i+1,j,k]!=0) and (bz[h,i,j,k]!=0)):
                        dbz_dx[h,i,j,k]=(bz[h,i+1,j,k]-bz[h,i,j,k])/(dx)

    dby_dz=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(k==0):
                        if((by[h,i,j,0]!=0) and (by[h,i,j,1]!=0)):
                            dbt_dy[h,i,j,0]=(by[h,i,j,1]-by[h,i,j,0])/dz
                    elif(k==nz-1):
                        if((by[h,i,j,nz-1]!=0) and (by[h,i,j,nz-2]!=0)):
                            dby_dz[h,i,j,nz-1]=(by[h,i,j,nz-1]-by[h,i,j,nz-2])/dz
                    elif((by[h,i,j,k-1]!=0) and (by[h,i,j,k+1]!=0)):
                        dby_dz[h,i,j,k]=(by[h,i,j,k+1]-by[h,i,j,k-1])/(2*dz)
                    elif((by[h,i,j,k-1]!=0) and (by[h,i,j,k]!=0)):
                        dby_dz[h,i,j,k]=(by[h,i,j,k]-by[h,i,j,k-1])/(dz)
                    elif((by[h,i,j,k+1]!=0) and (by[h,i,j,k]!=0)):
                        dby_dz[h,i,j,k]=(by[h,i,j,k+1]-by[h,i,j,k])/(dz)

    dbx_dz=np.zeros(temp.shape,dtype=np.float64)
    for h in range(nt):
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if(k==0):
                        if((bx[h,i,j,0]!=0) and (bx[h,i,j,1]!=0)):
                            dbx_dy[h,i,j,0]=(bx[h,i,j,1]-bx[h,i,j,0])/dz
                    elif(k==nz-1):
                        if((bx[h,i,j,nz-1]!=0) and (bx[h,i,j,nz-2]!=0)):
                            dbx_dz[h,i,j,nz-1]=(bx[h,i,j,nz-1]-bx[h,i,j,nz-2])/dz
                    elif((bx[h,i,j,k-1]!=0) and (bx[h,i,j,k+1]!=0)):
                        dbx_dz[h,i,j,k]=(bx[h,i,j,k+1]-bx[h,i,j,k-1])/(2*dz)
                    elif((bx[h,i,j,k-1]!=0) and (bx[h,i,j,k]!=0)):
                        dbx_dz[h,i,j,k]=(bx[h,i,j,k]-bx[h,i,j,k-1])/(dz)
                    elif((bx[h,i,j,k+1]!=0) and (bx[h,i,j,k]!=0)):
                        dbx_dz[h,i,j,k]=(bx[h,i,j,k+1]-bx[h,i,j,k])/(dz)
else:
    print("Error, method to compute derivative unknown...")
    sys.exit(2)


omega_tx=0.5*hbarc*(dbt_dx-dbx_dt)
omega_ty=0.5*hbarc*(dbt_dy-dby_dt)
omega_tz=0.5*hbarc*(dbt_dz-dbz_dt)

omega_yz=0.5*hbarc*(dby_dz-dbz_dy)
omega_zx=0.5*hbarc*(dbz_dx-dbx_dz)
omega_xy=0.5*hbarc*(dbx_dy-dby_dx)

with open(outputfile,"wb") as po:
     pickle.dump((tt,xx,yy,zz,vx,vy,vz,temp,omega_tx,omega_ty,omega_tz,omega_yz,omega_zx,omega_xy),po)
with open(outputfile+"_gradients","wb") as po:
     pickle.dump((tt,xx,yy,zz,vx,vy,vz,temp,bt,bx,by,bz,dbt_dx,dbt_dy,dbt_dz,dbx_dt,dby_dt,dbz_dt,dbx_dy,dbx_dz,dby_dx,dby_dz,dbz_dx,dbz_dy),po)
print("All done.")
