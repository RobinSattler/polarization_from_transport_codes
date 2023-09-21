# preprocess_thermodynamic_lattice_output_smash.py - version 0.1.1 - 04/01/2022

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import os.path
import glob
from timeit import default_timer as timer

sys.setrecursionlimit(10000)

#use binary (True) or ascii (False) input data
use_binary=True

#density tpye (it can be baryon or hadron)
density_type="hadron"
#density_type="baryon"

#if False it prints only error messages, if True it writes what it is doing at the moment and the intermediate results 
verbose=False

#output file version check
designed_lattice_version=1.0

if(verbose):
    init_start=timer()

#we parse the command line arguments
N_input_args=len(sys.argv)-1

if(N_input_args!=2):
   print ('Syntax: ./preprocess_thermodynamic_lattice_output_smash.py <data dir> <outputfile>')
   print ("Directory containing the files produced by SMASH with Thermodynamic Lattice Output")
   print ("outputfile is obviously the name of the output file with the results of the postprocessing")
   sys.exit(1)

#we get the name of input and output files
inputdir=sys.argv[1]
outputfile=sys.argv[2]

#we prepare lists of the input files
if (density_type == "hadron"):
    if(use_binary):
        net_bar_files=glob.glob(inputdir+'/hadron_j_QBS_*.bin')
        Tmunu_files=glob.glob(inputdir+'/hadron_tmn_landau_*.bin')
        vLandau_files=glob.glob(inputdir+'/hadron_v_landau_*.bin')
    else:
        net_bar_files=glob.glob(inputdir+'/hadron_j_QBS_*.dat')
        Tmunu_files=glob.glob(inputdir+'/hadron_tmn_landau_*.dat')
        vLandau_files=glob.glob(inputdir+'/hadron_v_landau_*.dat')
elif (density_type == "baryon"):
    if(use_binary):
        net_bar_files=glob.glob(inputdir+'/net_baryon_j_QBS_*.bin')
        Tmunu_files=glob.glob(inputdir+'/net_baryon_tmn_landau_*.bin')
        vLandau_files=glob.glob(inputdir+'/net_baryon_v_landau_*.bin')
    else:
        net_bar_files=glob.glob(inputdir+'/net_baryon_j_QBS_*.dat')
        Tmunu_files=glob.glob(inputdir+'/net_baryon_tmn_landau_*.dat')
        vLandau_files=glob.glob(inputdir+'/net_baryon_v_landau_*.dat')
else:
    print("Unknown density_type parameter (please, check the first lines of the script source code and fix it)")
    sys.exit(2)

#we sort the lists of input files
net_bar_files.sort()
Tmunu_files.sort()
vLandau_files.sort()

nf_bar=len(net_bar_files)
nf_tmn=len(Tmunu_files)
nf_vl=len(vLandau_files)

if((nf_bar != nf_tmn) or (nf_bar != nf_vl)):
    print("Sorry, but I can't continue.")
    print("I have found "+str(nf_bar)+" density current files, "+str(nf_tmn)+" Tmunu files, "+str(nf_vl)+" Landau velocity files")
    sys.exit(2)

#we check that all the input files have non zero length
for i in range(nf_bar):
    if((os.path.getsize(net_bar_files[i])==0) or (os.path.getsize(Tmunu_files[i])==0) or (os.path.getsize(vLandau_files[i])==0)):
        print("Because of a zero length file, I will not consider:")
        print(nf_bar.pop(i))
        print(nf_tmn.pop(i))
        print(nf_vl.pop(i))

#we update the length of the files
nf_bar=len(net_bar_files)
nf_tmn=len(Tmunu_files)
nf_vl=len(vLandau_files)

if(nf_bar*nf_tmn*nf_vl==0):
    print("Input files missing. I quit")
    sys.exit(2)
else:
    nf=nf_bar #we use a common variable for all file lengths

#dictionary containing information about the grid
lattice={}

#number of events
N_events=0

#correspondences between the indexes of the 1D Tmunu array and the energy momentum rank 2 tensor
iT00=0
iT01=1
iT02=2
iT03=3
iT10=1
iT11=4
iT12=5
iT13=6
iT20=2
iT21=5
iT22=7
iT23=8
iT30=3
iT31=6
iT32=8
iT33=9

#indexes corresponding to jB[0:3] in j_QBS
kb0=4
kb1=5
kb2=6
kb3=7

#empty list with the results
T_list=[]
jQBS_list=[]
v_list=[]

#empty list with the output times
tt=[]

#functions to read the header
def read_ascii_header(infile):
    lattice_version=float(infile.readline().split()[4])
    if(lattice_version != designed_lattice_version):
        print("Sorry, this code for analysis is designed for output version: "+str(designed_lattice_version)+", while you provided "+str(lattice_version)+"\n. I quit.")
        sys.exit(2)
    lattice_quantity=infile.readline().split()[1]
    lattice_dimensions=np.int32(infile.readline().split()[2:5])
    lattice_spacing=np.float64(infile.readline().split()[2:5])
    lattice_origin=np.float64(infile.readline().split()[2:5])
    return lattice_dimensions, lattice_spacing, lattice_origin

def read_binary_header(infile):
    lattice_version=np.fromfile(infile,dtype=np.float64,count=1)[0]
    if(lattice_version != designed_lattice_version):
        print("Sorry, this code for analysis is designed for output version: "+str(designed_lattice_version)+", while you provided "+str(lattice_version)+"\n. I quit.")
        sys.exit(2)
    lattice_quantity=np.fromfile(infile,dtype=np.int32,count=1)[0]
    lattice_dimensions=np.fromfile(infile,dtype=np.int32,count=3)
    lattice_spacing=np.fromfile(infile,dtype=np.float64,count=3)
    lattice_origin=np.fromfile(infile,dtype=np.float64,count=3)
    return lattice_dimensions, lattice_spacing, lattice_origin

def read_ascii_timestep_tmn(infile):
    #it returns True in case of EoF
    time_entry=infile.readline()
    if(time_entry==''):
        return True, None, None
    else:
        time=np.float64(time_entry)
    Tmunu=np.empty((10,nx,ny,nz),dtype=np.float64)
    for i in range(10):
        tmp_arr=np.loadtxt(infile,max_rows=ny*nz).flatten().reshape(nz,ny,nx).transpose()
        Tmunu[i,:,:,:]=tmp_arr.copy()
    return False, time, Tmunu

def read_binary_timestep_tmn(infile):
    #it returns True in case of EoF
    time_entry=np.fromfile(infile,dtype=np.float64,count=1)
    if(len(time_entry)==0):
        return True, None, None
    else:
        time=time_entry[0]
    Tmunu=np.empty((10,nx,ny,nz),dtype=np.float64)
    for i in range(10):
        tmp_arr=np.fromfile(infile,dtype=np.float64,count=nx*ny*nz).flatten().reshape(nz,ny,nx).transpose()
        Tmunu[i,:,:,:]=tmp_arr.copy()
    return False, time, Tmunu

def read_ascii_timestep_jqbs(infile):
    #we do not check the time entry
    time_entry=infile.readline()
    tmp_arr=np.loadtxt(infile,max_rows=nx*ny*nz).flatten().reshape(nz,ny,nx,12).transpose()
    return tmp_arr

def read_binary_timestep_jqbs(infile):
    #we do not check the time entry
    time_entry=np.fromfile(infile,dtype=np.float64,count=1)
    tmp_arr=np.fromfile(infile,dtype=np.float64,count=nx*ny*nz*12).flatten().reshape(nz,ny,nx,12).transpose()
    return tmp_arr

def read_ascii_timestep_vl(infile):
    #we do not check the time entry
    time_entry=infile.readline()
    tmp_arr=np.loadtxt(infile,max_rows=nx*ny*nz).flatten().reshape(nz,ny,nx,3).transpose()
    return tmp_arr

def read_binary_timestep_vl(infile):
    #we do not check the time entry
    time_entry=np.fromfile(infile,dtype=np.float64,count=1)
    tmp_arr=np.fromfile(infile,dtype=np.float64,count=nx*ny*nz*3).flatten().reshape(nz,ny,nx,3).transpose()
    return tmp_arr

for n_i in range(nf):
    i_tmn=Tmunu_files[n_i]
    i_jqbs=net_bar_files[n_i]
    i_vl=vLandau_files[n_i]
    if(verbose):
        print("Opening "+i_tmn)
        start_time = timer()
    if(use_binary):
        if(i_tmn[-3:]==".gz"):
            fp_tmn=gzip.open(i_tmn,"rb")
        else:
            fp_tmn=open(i_tmn,"rb")
        lattice_dimensions, lattice_spacing, lattice_origin = read_binary_header(fp_tmn)
    else:
        if(i_tmn[-3:]==".gz"):
            fp_tmn=gzip.open(i_tmn,"r")
        else:
            fp_tmn=open(i_tmn,"r")
        lattice_dimensions, lattice_spacing, lattice_origin = read_ascii_header(fp_tmn)
    if(verbose):
        print("Opening "+i_jqbs)
        start_time = timer()
    if(use_binary):
        if(i_jqbs[-3:]==".gz"):
            fp_jqbs=gzip.open(i_jqbs,"rb")
        else:
            fp_jqbs=open(i_jqbs,"rb")
        lattice_dimensions_jqbs, lattice_spacing_jqbs, lattice_origin_jqbs = read_binary_header(fp_jqbs)
    else:
        if(i_jqbs[-3:]==".gz"):
            fp_jqbs=gzip.open(i_jqbs,"r")
        else:
            fp_jqbs=open(i_jqbs,"r")
        lattice_dimensions_jqbs, lattice_spacing_jqbs, lattice_origin_jqbs = read_ascii_header(fp_jqbs)
    if(verbose):
        print("Opening "+i_vl)
        start_time = timer()
    if(use_binary):
        if(i_vl[-3:]==".gz"):
            fp_vl=gzip.open(i_vl,"rb")
        else:
            fp_vl=open(i_vl,"rb")
        lattice_dimensions_vl, lattice_spacing_vl, lattice_origin_vl = read_binary_header(fp_vl)
    else:
        if(i_vl[-3:]==".gz"):
            fp_vl=gzip.open(i_vl,"r")
        else:
            fp_vl=open(i_vl,"r")
        lattice_dimensions_vl, lattice_spacing_vl, lattice_origin_vl = read_ascii_header(fp_vl)

    #we skip the check that the grid data for the various quantities are compatible with each other 

    if(n_i==0): #initially we need to acquire some information
        lattice["dimensions"]=lattice_dimensions
        lattice["spacing"]=lattice_spacing
        lattice["origin"]=lattice_origin
        nx,ny,nz=lattice_dimensions[:]
    else:
        if(not np.array_equal(lattice["dimensions"],lattice_dimensions)):
            print("Error in file "+i_tmn+": different lattice dimensions. Until now: "+str(lattice["dimensions"])+", this time: "+str(lattice_dimensions)+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice["spacing"],lattice_spacing)):
            print("Error in file "+i_tmn+": different lattice spacing. Until now: "+str(lattice["spacing"])+", this time: "+str(lattice_spacing)+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice["origin"],lattice_origin)):
            print("Error in file "+i_tmn+": different lattice origin. Until now: "+str(lattice["origin"])+", this time: "+str(lattice_origin)+".\nI quit.")
            sys.exit(2)

    index=0

    while(True):
        if(use_binary):
            eof,time,Tmunu=read_binary_timestep_tmn(fp_tmn)
            if(eof):
                break
            j_QBS=read_binary_timestep_jqbs(fp_jqbs) #first index: j component then x, y, z
            vl=read_binary_timestep_vl(fp_vl) #first index: v component then x, y, z
        else:
            eof,time,Tmunu=read_ascii_timestep_tmn(fp_tmn)
            if(eof):
                break
            j_QBS=read_ascii_timestep_jqbs(fp_jqbs) #first index: j component then x, y, z
            vl=read_ascii_timestep_vl(fp_vl) #first index: v component then x, y, z
        if(n_i==0):
            tt.append(time)
        else:
            if(time!=tt[index]):
                print("Error when reading file "+i+":")
                print("At step "+str(index)+" time "+str(tt[index])+" was expected, while "+str(time)+" was found. I quit.\n")
                sys.exit(2)

        if(n_i==0):
            T_list.append(Tmunu)
            jQBS_list.append(j_QBS)
            v_list.append(vl)
        else:
            T_list[index]=T_list[index]+Tmunu
            jQBS_list[index]=jQBS_list[index]+j_QBS
            v_list[index]=v_list[index]+vl

        index=index+1

        if(verbose):
            print("Done timestep: "+str(index)+", simulation time: "+str(time))

    #we check that we counted correctly the timesteps:
    if(len(tt)!=index):
        print("Error, I counted "+str(index)+" timesteps, but I have "+str(len(tt))+" entries in the list of timesteps...\nI quit.")
        sys.exit(2)
    if(n_i==0):
        nt=index

    fp_tmn.close()
    fp_jqbs.close()
    fp_vl.close()
    if(verbose):
        end_time = timer()
        print(i_tmn+", "+i_jqbs+", "+i_vl+" read in "+tf.format(end_time-start_time)+" seconds")
N_events=n_i+1

# we transform the lists of 3D arrays into 4D arrays
# we do it one by one to save memory
T_arr=np.zeros((nt,10,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    T_arr[i,:,:,:]=T_list[i][:,:,:,:]
T_list=None
del T_list
jQBS_arr=np.zeros((nt,12,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    jQBS_arr[i,:,:,:]=jQBS_list[i][:,:,:,:]
jQBS_list=None
del jQBS_list
v_arr=np.zeros((nt,3,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    v_arr[i,:,:,:]=v_list[i][:,:,:,:]
v_list=None
del v_list


if(verbose):
    print("Writing the final results in "+outputfile)
    start_time = timer()

with open(outputfile,"wb") as po:
      pickle.dump((lattice,tt,N_events,T_arr,jQBS_arr,v_arr),po)

if(verbose):
    end_time = timer()
    print("Done in "+tf.format(end_time-start_time)+" seconds")
    print("All done in "+tf.format(end_time-init_start)+" seconds")
