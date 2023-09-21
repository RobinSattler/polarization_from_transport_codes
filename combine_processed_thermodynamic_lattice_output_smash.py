# combine_processed_thermodynamic_lattice_output_smash.py - version 0.1.2- 04/01/2022

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import gzip
from timeit import default_timer as timer

#if False it prints only error messages, if True it writes what it is doing at the moment and the intermediate results 
verbose=True

if(verbose):
    init_start=timer()

#we parse the command line arguments
N_input_args=len(sys.argv)-1

if(N_input_args<2):
   print ('Syntax: ./combine_processed_thermodynamic_lattice_output_smash.py <file data 1> [data 2] ... <outputfile>')
   print ("file data 1,2,3...N are the pickled files produced by preprocess_thermodynamic_lattice_output_smash.py")
   print ("outputfile is obviously the name of the output file with the results of the postprocessing")
   sys.exit(1)

#we get the name of input and output files
inputfiles=sys.argv[1:N_input_args]
n_input=len(inputfiles)
outputfile=sys.argv[N_input_args]

for n_i, infile in enumerate(inputfiles):
    if(infile[-3:]==".gz"):
       if(verbose):
           print("Opening gzipped file "+infile)
       pi=gzip.open(infile,"rb")
    else:
        if(verbose):
             print("Opening file "+infile)
        pi=open(infile,"rb")

    indata=pickle.load(pi)
    pi.close()
    lattice,tt,N_file_events,T_arr,jQBS_arr,v_arr = indata[:]

    if (verbose):
        print("Adding "+str(N_file_events)+" events")
    if(n_i==0):
        lattice_ref=lattice
        tt_ref=tt
        nt=len(tt)
        Tmunu=T_arr.copy()
        jQBS=jQBS_arr.copy()
        v=v_arr.copy()
        N_events=N_file_events
    else:
        if(not np.array_equal(lattice_ref["dimensions"],lattice["dimensions"])):
            print("Error in file "+infile+": different lattice dimensions. Until now: "+str(lattice_ref["dimensions"])+", this time: "+str(lattice["dimensions"])+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice_ref["spacing"],lattice["spacing"])):
            print("Error in file "+infile+": different lattice spacing. Until now: "+str(lattice_ref["spacing"])+", this time: "+str(lattice["spacing"])+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice_ref["origin"],lattice["origin"])):
            print("Error in file "+infile+": different lattice origin. Until now: "+str(lattice_ref["origin"])+", this time: "+str(lattice["origin"])+".\nI quit.")
            sys.exit(2)
            sys.exit(2)
        N_events=N_events+N_file_events
        Tmunu=Tmunu+T_arr
        jQBS=jQBS+jQBS_arr
        v=v+v_arr


if (verbose):
    print ("Writing output file "+outputfile+" based on "+str(N_events)+" events")
with open(outputfile,"wb") as po:
    pickle.dump((lattice,tt,N_events,Tmunu,jQBS,v),po)
