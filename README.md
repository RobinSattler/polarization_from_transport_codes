## Information about the scripts contained in this repository
### The exact syntax for each script can be obtained by invoking it without any argument
### In the following explanations
### nt is the total number of timesteps
### nx, ny and nz are the number of cells in the x, y and z direction, respectively
### tmn is the energy momentum tensor in the Landau frame
### v is the fluid velocity (i.e. the velocity of the Landau frame in the computational frame)
### jBQS are the net baryon, electric charge and strangeness currents, respectively

* A - preprocess_thermodynamic_lattice_output_smash.py
  
  It takes as argument a directory in which the thermodynamic lattice data files are stored and
  produces a python pickle archive containing the following data: the structure of the lattice,
  and array with the values of the timesteps, the number of events, tmp as a numpy array with shape
  nt,10,nx,ny,nz, jQBS as a numpy array with shape nt,12,nx,ny,nz, v as a numpy array with shape
  nt,3,nx,ny,nz. See the code for more details about the internal representation.
  By editing the first line of the script it is possible to choose as input format either the ascii
  or the binary thermodynamic lattice SMASH output. It is also possible to choose the density type
  (hadron or baryon).

* B - combine_processed_thermodynamic_lattice_output_smash.py

  It combines multiple outputs of the script A into a single file with the same structure.

* C - compute_vorticity_cg_data.py

  It takes the output of UrQMD data processed with another script not in this repository
  and computes the thermal vorticity.
  
* D - compute_vorticity_from_th_latt_output_smash.py
  
  It takes the output of A or B and it computes the thermal vorticity.
  The script needs also a tabulated EoS, which can be either the UrQMD or the SMASH HG EoS.
  The data about the EoS are hardcoded at the beginning of the script.
  The EoS are not provided in this repository.
  The script produces two pickle archive files: one with the vorticity components and one with
  just the partial derivatives.

* E - make_vorticity_plots_smash.py

  It plots the data of the vorticity components produced by the script D.

* F - make_vort_deriv_plots_for_dbg_smash.py

  It plots the data of the partial derivatives produced by the script D.

* G - compute_mean_spin_smash.py

  It takes in input the output of D and a list of position and momenta of hadrons of the same species and
  produces in output a list of the same hadrons with the components of their mean spin.
  The list with the coordinates and the momenta of the hadrons at the moment of their chemical and kinetic
  freezeout is created by a script not included in this repository from SMASH data.
  By editing a few lines in the first part of the script it is possible to choose which kind of freezeout
  data must be used and to introduce transverse momentum and rapidity cuts.

* H - compute_mean_spin_smash_Oscar_GM_files.py

  It is similar to the script G, but the of list of position and momenta of hadrons that takes in input
  must be produced by a script written by Oscar Garcia Montero.

* I - compute_mean_spin_urqmd.py
 
  It takes in input the output of C and a list of position and momenta of hadrons and produces in output  
  lists of the same hadrons (one for each hadron type) with the components of their mean spin.
  The list with the coordinates and the momenta of the hadrons at the moment of their chemical or kinetic 
  freezeout is created by a script not included in this repository from UrQMD data.
  By editing a few lines in the first part of the script it is possible to choose which hadrons species
  to select and to introduce transverse momentum and rapidity cuts.

* L - compute_dN_dt_from_urqmd_hadron_data.py
 
  It takes the same hadron list used by script I and it produces in output text files with their
  dN/dt (t) distribution. Currently the script works with Lambda and Sigma baryons only, albeit it can be
  easily extended.

* M - compute_dN_dt_dN_dxdz_from_urqmd_hadron_data.py

  It takes in input the output of the script I and it produces in output text files with time and spatial
  distribution (on the zx plane) of the hadrons at the various timesteps.

* N - compute_dN_vs_time_and_space.py

  It takes in input the output of the script G and it produces in output text files with the time and space
  distributions of the hadron species in the input file.
