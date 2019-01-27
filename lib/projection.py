"""
  Projects out translation, external rotation, and hindered rotor motions from the Hessian
"""

import sys
import numpy as np

import molecule

def project_Hessian(masses, xyz, H):
    """ Project the motions out of the Hessian. """

    # Calculate the vectors associated with translation, external rotation, and hindered internal rotation 
    TransMat   = build_translation_vectors(masses)
    ExtRotMat  = build_ext_rotation_vectors(masses, xyz)
    # HindRotMat = build_hind_rotation_vectors()
    NHindRot = 0
    Natoms = len(mass)

    # Compound the motional vectors into a vector D
    D = np.column_stack( (TransMat, ExtRotMat) )

    # print(TransMat)
    print(ExtRotMat)
    # print(HindRotMat)
    # print(D)

    # Orthogonalize D matrix
    # Dortho = orthogonalize_D_matrix(D)
  
    # Build the Projection Matrix (P)
    Ds = np.dot(D, D.T)
    #print(Ds)
    I = np.identity(3*Natoms)
    P = I - Ds

    # Project out the modes out of the Hessian Matrix 
    # Hproj = np.dot(P.T * H * P)

    # Diagonalize the Projected Hessian
    # lam, L = np.linalg.eig(Hproj)
    lam = 1.0

    return lam


def build_translation_vectors(masses):
    """ Builds the vectors that describe the translation of the molecule.
    """    

    # Calculate number of atoms
    Natoms = len(masses)

    # Build the vectors for translation
    Transx = np.zeros(3*Natoms) 
    Transy = np.zeros(3*Natoms) 
    Transz = np.zeros(3*Natoms) 
    for i in range(Natoms):
        Transx[3*i:3*i+3] = tuple(np.sqrt(masses[i]) * elem for elem in [1.0, 0.0, 0.0])
        Transy[3*i:3*i+3] = tuple(np.sqrt(masses[i]) * elem for elem in [0.0, 1.0, 0.0])
        Transz[3*i:3*i+3] = tuple(np.sqrt(masses[i]) * elem for elem in [0.0, 0.0, 1.0])

    # Put vectors into a matrix
    Trans = np.column_stack( (Transx, Transy, Transz) )

    return Trans


def build_ext_rotation_vectors(masses, xyz):
    """ Builds the vectors that describe the external rotation of the molecule. 
    """    

    # Calculate number of atoms
    Natoms = len(masses)
 
    # Calculate the moment-of-inertia and center-of-mass
    Rcom = molecule.calc_center_of_mass(masses, xyz)
    Ip, Xp = molecule.calc_moment_of_inertia(masses, xyz)

    # Obtain vectors for the projection of ...
    Px = np.dot(Rcom, X[0])
    Py = np.dot(Rcom, X[1])
    Pz = np.dot(Rcom, X[2])

    RotExt = np.zeros((3*Natoms,3),dtype=float)
    for i in range(3*Natoms):
        for j in range(3):
            RotExt[i][j] = (Py[i]*Xp[j][2] - Pz[i]*Xp[j][1]) / np.sqrt(masses[i])
            RotExt[i][j] = (Pz[i]*Xp[j][0] - Px[i]*Xp[j][2]) / np.sqrt(masses[i])
            RotExt[i][j] = (Px[i]*Xp[j][1] - Py[i]*Xp[j][0]) / np.sqrt(masses[i])

    return RotExt


def build_hind_rotor_vectors(xyz):
    """ Builds the vectors that describe the rotaiton of a hindered rotor. 
    """    

    # Initialize list to hold all of the Rotor vectors
    Nrotors = 0
    Rotors = []

    for i in range(Nrotors):

        # Numbers of the pivot atoms
        pivotA = 
        pivotB =   
  
    # Numbers denoting atoms attached to pivot atoms A and B 
    attA = []
    attB = [] 
  
    for j in range(len(attA)):
      np.cross( (xyz[j]-xyz[pivotA]), (xyz[pivotA]-xyz[pivotB]) ) / np.norm(xyz[pivotA]-xyz[pivotB]) 
    for j in range(len(attB)):
      np.cross( (xyz[j]-xyz[pivotB]), (xyz[pivotB]-xyz[pivotA]) ) / np.norm(xyz[pivotB]-xyz[pivotA]) 

  return None 


def ortho_D_matrix():
  """ Orthogonalize the matrix describing the projected motions. 
  """
  return None


#if __name__ == 'main':
REAC, TS, temperature = molecule.parse_input(sys.argv[1])
print(project_Hessian(TS['masses,'], TS['xyz'], TS['hessian']))

##### END PROGRAM #####


