'''
  Projects out translation, external rotation, and hindered rotor motions from the Hessian
'''

import sys
import numpy as np
import math
import molecule
from constants import atom_mass

def project_Hessian(geom):
  ''' Project the motions out of the Hessian. '''

  # Get the masses
  atom_symbols = []
  for coord in geom:
    atom_symbols.append( coord[0] )
  mass = np.asarray( [ atom_mass[symbol] for symbol in atom_symbols ] ) 

  # Calculate the vectors associated with translation, external rotation, and hindered internal rotation 
  TransMat   = build_translation_vectors(mass)
  ExtRotMat  = build_ext_rotation_vectors(mass, geom)
  #HindRotMat = build_hind_rotation_vectors()
  NHindRot = 0
  Natoms = len(mass)

  # Compound the motional vectors into a vector D
  D = np.column_stack( (TransMat, ExtRotMat) )

  #print(TransMat)
  print(ExtRotMat)
  #print(HindRotMat)
  #print(D)

  # Orthogonalize D matrix
  #Dortho = orthogonalize_D_matrix(D)
  
  # Build the Projection Matrix (P)
  Ds = np.dot(D, D.T)
  #print(Ds)
  I = np.identity(3*Natoms)
  P = I - Ds

  # Project out the modes out of the Hessian Matrix 
  #Hproj = np.dot(P.T * H * P

  # Diagonalize the Projected Hessian
  #lam, L = np.linalg.eig(Hproj)
  lam = 1.0

  return lam


def build_translation_vectors(masses):
  ''' Builds the vectors that describe the translation of the molecule.
  '''    

  # Calculate number of atoms
  Natoms = len(masses)

  # Build the vectors for translation
  Transx, Transy, Transz = np.zeros(3*Natoms), np.zeros(3*Natoms), np.zeros(3*Natoms) 
  for i in range(Natoms):
    Transx[ 3*i : 3*i+3 ] = tuple( np.sqrt(masses[i]) * elem for elem in [ 1.0, 0.0, 0.0 ] )
    Transy[ 3*i : 3*i+3 ] = tuple( np.sqrt(masses[i]) * elem for elem in [ 0.0, 1.0, 0.0 ] )
    Transz[ 3*i : 3*i+3 ] = tuple( np.sqrt(masses[i]) * elem for elem in [ 0.0, 0.0, 1.0 ] )

  # Put vectors into a matrix
  Trans = np.column_stack( (Transx, Transy, Transz) )

  return Trans


def build_ext_rotation_vectors(masses, geom):
  ''' Builds the vectors that describe the external rotation of the molecule. 
  '''    

  # Calculate number of atoms
  Natoms = len(masses)
 
  # Calculate the moment-of-inertia and center-of-mass
  mom_inert, X, Rcom = molecule.calc_moment_of_inertia(geom)

  Px = np.dot( Rcom, X[0] )
  Py = np.dot( Rcom, X[1] )
  Pz = np.dot( Rcom, X[2] )

  RotExt = np.zeros((3*Natoms,3),dtype=float)
  for i in range(3*Natoms):
    for j in range(3):
      RotExt[i][j] = ( Py[i]*X[j][2] - Pz[i]*X[j][1] ) / np.sqrt(masses[i])
      RotExt[i][j] = ( Pz[i]*X[j][0] - Px[i]*X[j][2] ) / np.sqrt(masses[i])
      RotExt[i][j] = ( Px[i]*X[j][1] - Py[i]*X[j][0] ) / np.sqrt(masses[i])

  return RotExt


def build_hind_rotor_vectors():
  ''' Builds the vectors that describe the rotaiton of a hindered rotor. 
  '''    
  
  
  
  return None 


def ortho_D_matrix():
  ''' Orthogonalize the matrix describing the projected motions. 
  '''
  return None


#if __name__ == 'main':
REAC, TS, temperature = molecule.parse_input(sys.argv[1])
print(project_Hessian(REAC['geom']))

##### END PROGRAM #####


