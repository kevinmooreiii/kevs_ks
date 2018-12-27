'''
  Parses the input file to get user-specified information about the molecular species and calculates other structural information.
'''

import numpy as np
from constants import atom_mass


def parse_input(input_file_name):
  ''' Read the input file and get the info 
  '''  

  # Open the file and read the lines
  with open(input_file_name, 'r') as infile:
    lines = infile.readlines()  

  # Obtain the range of temperatures of energies to for which to calculate rate constants 
  k_var = []
  for i in range(len(lines)):
    if 'temperature =' in linesor 'energies =' in lines[:
      tmp = line.strip().split()
      k_var = range(tmp[2], tmp[3], tmp[4])
      break

  # Determine the the line numbers in the input file where the reactant, TS, and product info stored
  reac1_start, reac1_end = 0, 0
  #reac2_start, reac2_end = 0, 0
  ts_start, ts_end = 0, 0
  #prod1_start, prod1_end = 0, 0
  #prod2_start, prod2_end = 0, 0
  for i in range(len(lines)):
    if '$REAC1' in lines[i]:
      reac1_start = i
      for j in range(reac1_start, len(lines)):
        if '$END' in lines[j]:
          reac1_end = j
    #if '$REAC2' in lines[i]:
    #  reac2_start = i
    #  for j in range(reac2_start, len(lines)):
    #    if '$END' in lines[j]:
    #      reac2_end = j
    if '$TS' in lines[i]:
      ts_start = i
      for j in range(ts_start, len(lines)):
        if '$END' in lines[j]:
          ts_end = j
    #if '$PROD1' in lines[i]:
    #  prod1_start = i
    #  for j in range(prod1_start, len(lines)):
    #    if '$END' in lines[j]:
    #      prod1_end = j
    #if '$PROD2' in lines[i]:
    #  prod2_start = i
    #  for j in range(prod2_start, len(lines)):
    #    if '$END' in lines[j]:
    #      prod2_end = j

  # Grab the reactant information
  reac1_moltype = lines[ts_start].strip().split()[1].upper() 
  reac1_masses, reac1_coords, reac1_energy, reac1_freqs, reac1_sym, reac1_eleclvls, reac_tunn = [], '', [], '', [], []
  for i in range(reac_start, reac_end):

    if 'Geometry' in lines[i]:
      geom_start = i + 1
      geom_end = geom_start + int(lines[i].strip().split()[1])
      reac1_masses = [ atom_mass[lines[j].strip.split()[0]] for j in range(geom_start, geom_end) ]     
      reac1_coords = [ lines[j].strip().split()[1:] for j in range(geom_start, geom_end) ]     

    if 'Energy' in lines[i]:
      reac1_energy = float(lines[i+1].strip())

    if 'Frequencies' in lines[i]:
      freq_start = i + 1
      freq_end = freq_start + int(lines[i].strip().split()[1])
      reac1_freqs = [ float(lines[j].strip.split()) for j in range(geom_start, geom_end) ]     

    if 'Symmetry' in lines[i]:
      reac1_sym = lines[i+1].strip().upper()

    if 'ElecLevels' in lines[i]:
      eleclvls_start = i + 1
      eleclvls_end = eleclvls_start + int(lines[i].strip().split()[1])
      reac1_eleclvls = [ lines[j].strip().split()[1:] for j in range(geom_start, geom_end) ]     
    
  # Grab the transition state information
  ts_moltype = lines[ts_start].strip().split()[1].upper() 
  ts_masses, ts_coords, ts_energy, ts_freqs, ts_sym, ts_eleclvls, ts_tunn = [], '', [], '', [], []
  for i in range(ts_start, ts_end):

    if 'Geometry' in lines[i]:
      geom_start = i + 1
      geom_end = geom_start + int(lines[i].strip().split()[1])
      ts_masses = [ atom_mass[lines[j].strip.split()[0]] for j in range(geom_start, geom_end) ]     
      ts_coords = [ lines[j].strip().split()[1:] for j in range(geom_start, geom_end) ]     

    if 'Energy' in lines[i]:
      reac_energy = float(lines[i+1].strip())

    if 'Frequencies' in lines[i]:
      freq_start = i + 1
      freq_end = freq_start + int(lines[i].strip().split()[1])
      ts_freqs = [ float(lines[j].strip.split()) for j in range(geom_start, geom_end) ]     

    if 'Symmetry' in lines[i]:
      ts_sym = lines[i+1].strip().upper()

    if 'ElecLevels' in lines[i]:
      eleclvls_start = i + 1
      eleclvls_end = eleclvls_start + int(lines[i].strip().split()[1])
      ts_eleclvls = [ lines[j].strip().split()[1:] for j in range(geom_start, geom_end) ]     
    
    if 'Tunneling' in lines[i]:
      ts_tunn =  lines[i+1].strip().split() 

  # Store the information in dictionaries  
  REAC1 = {
    'moltype'    : reac1_moltype,
    'masses'     : reac1_masses,
    'xyz'        : reac1_xyz,
    'energy'     : reac1_energy,
    'freqs'      : reac1_freqs,
    'sym'        : reac1_sym,
    'eleclevels' : reac1_eleclevels,
  }

  TS = {
    'moltype'    : ts_moltype,
    'masses'     : ts_masses,
    'xyz'        : ts_xyz,
    'energy'     : ts_energy,
    'freqs'      : ts_freqs,
    'sym'        : ts_sym,
    'eleclevels' : ts_eleclevels,
    'tunneling'  : ts_tunn 
  }

  return REAC, TS, k_var


def det_rotation_number(sym):
  ''' Determine the rotational number using the molecular symmetry. '''   
  
  sigma = constants.rot_sym_num[sym]

  return sigma


def calc_moment_of_inertia(masses, xyz):
  ''' Find the principal moments of intertia using the molecular xyz coordinates. '''

  # Determine the coordinates of the center-of-mass
  xm, ym, zm, mass_total = 0.0, 0.0, 0.0, 0.0
  for i in range(len(masses)):
    xm += masses[i] * xyz[i,0]    
    ym += masses[i] * xyz[i,1]    
    zm += masses[i] * xyz[i,2]    
    mass_total += masses[i]
  
  R_com = np.array( [ (xm / mass_total), (ym / mass_total), (zm / mass_total) ] )
  
  # Translate the coordinates to the center-of-mass
  R_trans = xyz - R_com

  # Build the moment-of-inertia tensor
  Ixx, Iyy, Izz, Ixy, Iyx, Ixz, Izx, Iyz, Izy = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  
  for i in range( np.shape(xyz_trans)[0] ):
    Ixx += masses[i] * ( R_trans[i][1]**2 + R_trans[i][2]**2 )  
    Iyy += masses[i] * ( R_trans[i][0]**2 + R_trans[i][2]**2 )  
    Izz += masses[i] * ( R_trans[i][0]**2 + R_trans[i][1]**2 )  
    Ixy += -1.0 * masses[i] * ( R_trans[i][0] * R_trans[i][1] )  
    Iyx += -1.0 * masses[i] * ( R_trans[i][0] * R_trans[i][1] )  
    Ixz += -1.0 * masses[i] * ( R_trans[i][0] * R_trans[i][2] )  
    Izx += -1.0 * masses[i] * ( R_trans[i][0] * R_trans[i][2] )  
    Iyz += -1.0 * masses[i] * ( R_trans[i][1] * R_trans[i][2] )  
    Izy += -1.0 * masses[i] * ( R_trans[i][1] * R_trans[i][2] )  
  
  I = np.array( [ [ Ixx, Ixy, Ixz],
                  [ Iyx, Iyy, Iyz],
                  [ Izx, Izy, Izz] ] )

  # Diagonalize the moment-of-inertia tensor(I) to get 
  # moments-of-inertia () and the principal axes of rotation (X)
  Ip, Xp = np.linalg.eig(I) 

  return Ip, Xp, R_trans

