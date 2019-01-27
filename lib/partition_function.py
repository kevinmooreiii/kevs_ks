'''
Computes the vibrational, rotational, translational, and electronic partition functions for a given molecule
'''

import numpy as np

import molecule
import projection
from constants import h, kB, atom_mass
from constants import cm_to_Hz, amu_to_kg, ang_to_m


def calc_Q_total(T, moltype, masses, geom, energy, freqs, sym, eleclevels):
  ''' Calculates the total partition function using the components. 
  '''

  # Calculate the components of the total partition function Q 
  q_vib   = calc_q_vib(freqs, T)  
  q_rot   = calc_q_rot(masses, xyz, sym, moltype, T)  
  q_trans = calc_q_trans(geom, T)  
  q_elec  = calc_q_elec(eleclevels, T)  

  # Assume the components of the total partition function are seperable 
  Q_total = q_vib * q_rot * q_trans * q_elec

  return Q_total


def calc_q_vib(freqs, T):
  ''' Calculates the vibrational partition function. 
      Harmonic Oscillator approximation used. 
  '''
    
  # Simply use the harmonic frequencies to compute partition function in absense of hindered rotors
  if Nrotors == 0:

    q_vib = 1.0
    for v in freqs:
      v = v * cm_to_Hz
      q_vib *= ( ( np.exp( (-h*v) / (2.0*kB*T) ) ) / ( 1.0 - np.exp( (-h*v) / (kB*T) ) ) )
 

  # Project the Hessian to obtain projected harmonic freqeuncies and hindered rotors to calculate the partition functions 
  else:

    # Take the provided Hessian and project out the translations, external rotations, and hindered rotors
    # Diagonalize the Hessian, obtain the harmonic frequencies, compute the partition function for these vibrations
    proj_freqs = calc_projected_freqs(hessian, rotor_stuff)
    
    # Calculate the partition function for the projected vibrational frequencies
    q_proj_freqs = 1.0
    for v in proj_freqs:
      v = v * cm_to_Hz
      q_vib_projfreqs *= ( ( np.exp( (-h*v) / (2.0*kB*T) ) ) / ( 1.0 - np.exp( (-h*v) / (kB*T) ) ) )
    
    # Calculate the partition function for the hindered rotors
    q_hind_rots = 1.0
    for n in Nrotors:
      q_vib_hindrots *= calc_q_hind_rots
    
    # Calculate the total vibrational partition function   
    q_vib = q_vib_projfreqs * q_vib_hindrots


  return q_vib


def calc_q_rot(masses, xyz, sym, moltype, T):
  ''' Calculates the rotational partition function. 
      Rigid-Rotor approximation used.
  '''

  if moltype == 'ATOM': 
    
    # Set partition function equal to 1 for an atom
    q_rot = 1.0
 
  elif moltype == 'LINEAR' or moltype == 'NONLINEAR':

    # Determine the rotational symmetry number
    sigma = molecule.det_rotation_number(sym)

    # Calculate principal moments-of-intertia
    Ip, Xp = molecule.calc_moment_of_inertia(masses, xyz)
    
    # Compute q_rot for linear or polyatomics
    if moltype == 'LINEAR': 
      q_rot = ( 1.0 / sigma ) * ( (8.0*(np.pi**2)*I*kB*T) / h**2 )
    elif moltype == 'NONLINEAR':
      q_rot = np.sqrt(np.pi) / sigma
      for I in Ip:
        I = I * amu_to_kg * ang_to_m * ang_to_m 
        q_rot *= np.sqrt( (8.0*(np.pi**2)*I*kB*T) / h**2 )

  return q_rot


def calc_q_trans(masses, T):
  ''' Calculates the translational partition function. 
      Ideal Gass approximation used.
  '''

  # Sum the masses of the atoms to get the total mass of the molecule
  Natoms = len(masses)
  for i in range(Natoms):
    Mass += masses[i] * amu_to_kg 
  
  # Compute the partition function 
  q_trans = ( (2.0*np.pi*Mass*kB*T) / (h**2) )**(3.0/2.0) * ( (kB*T) / 101325.0 ) 
  
  return q_trans


def calc_q_elec(eleclevels, T):
  ''' Calculates the electronic partition function. 
  '''

  # Loop over the energy levels
  q_elec = 0.0
  for level in eleclevels:
    q_elec += level[0] * np.exp( (-level[1]) / (kB*T) )

  return q_elec


def calc_q_hindrot():

  # Integrate over the hindered rotor potential to get the partition function
  q_hindrot = integral()

  return q_hindrot


