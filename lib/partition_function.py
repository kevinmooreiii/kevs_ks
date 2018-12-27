'''
Computes the vibrational, rotational, translational, and electronic partition functions for a given molecule
'''

import math
import numpy as np
import molecule
#import constants 
from constants import h, kB, atom_mass
from constants import cm_to_Hz, amu_to_kg, ang_to_m


def calc_Q_total(T, moltype, geom, energy, freqs, sym, eleclevels):
  ''' Calculates the total partition function using the components'''

  # Calculate the components of the total partition function Q 
  q_vib   = calc_q_vib(freqs, T)  
  q_rot   = calc_q_rot(geom, sym, moltype, T)  
  q_trans = calc_q_trans(geom, T)  
  q_elec  = calc_q_elec(eleclevels, T)  

  # Assume the components of the total partition function are seperable 
  Q_total = q_vib * q_rot * q_trans * q_elec

  return Q_total


def calc_q_vib(freqs, T):
  ''' Calculates the vibrational partition function. '''

  q_vib = 1.0
  for v in freqs:
    v = v * cm_to_Hz
    q_vib *= ( ( np.exp( (-h*v) / (2.0*kB*T) ) ) / ( 1.0 - np.exp( (-h*v) / (kB*T) ) ) )

  return q_vib


def calc_q_rot(geom, sym, moltype, T):
  ''' Calculates the rotational partition function. '''

  if moltype == 'ATOM': 
    
    # Set partition function equal to 1 for an atom
    q_rot = 1.0
 
  else:

    # Determine the rotational symmetry number
    sigma = molecule.det_rotation_number(sym)

    # Calculate principal moments-of-intertia
    mom_inert, vec = molecule.calc_moment_of_inertia(geom)
    
    # Compute q_rot for linear or polyatomics
    if moltype == 'LINEAR': 
      q_rot = ( 1.0 / sigma ) * ( (8.0*(np.pi**2)*I*kB*T) / h**2 )
    elif moltype == 'NONLINEAR':
      q_rot = np.sqrt(np.pi) / sigma
      for I in mom_inert:
        I = I * amu_to_kg * ang_to_m * ang_to_m 
        q_rot *= np.sqrt( (8.0*(np.pi**2)*I*kB*T) / h**2 )

  return q_rot


def calc_q_trans(geom, T):
  ''' Calculates the translational partition function. '''

  # Get the mass
  Natoms = len(geom)
  Mass = 0.0
  for i in range(Natoms):
    Mass += atom_mass[geom[i][0]] * amu_to_kg 
  
  q_trans = ( (2.0*np.pi*Mass*kB*T) / (h**2) )**(3.0/2.0) * ( (kB*T) / 101325.0 ) 
  
  return q_trans


def calc_q_elec(eleclevels, T):
  ''' Calculates the electronic partition function. '''

  # Initialize electronic partition function to zero
  q_elec = 0.0

  # Loop over the energy levels
  for level in eleclevels:
    q_elec += float(level[0]) * np.exp( (-float(level[1])) / (kB*T) )

  return q_elec


