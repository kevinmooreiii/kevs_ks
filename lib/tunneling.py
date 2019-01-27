'''
  Calculates coefficients to account for quantum tunneling effects to the rate constants.
  Tunneling Models Supported: NoTunnel, Eckart, and Wigner 
'''

from constants import kB, h


def calc_tunn_coeff(E_reac, E_prod, model, E_TS, vTS, T):
  ''' Calls the correct tunneling correction based on the tunneling model.
  '''
  
  if model == 'NOTUNNEL':
    kappa_T = 1.0
  elif model == 'ECKART':
    kappa_T = calc_eckart_tunn_coeff(E_reac, E_prod, E_TS, vTS, T)
  elif model == 'WIGNER':
    kappa_T = calc_wigner_tunn_coeff(vTS, T)

  return kappa_T


def calc_eckart_tunn_coeff(E_reac, E_prod, E_TS, vTS, T):
  ''' Calculates kappa(T) using the Eckart potential model. 
  '''

  # Set the potential values so that the min E is set to the side with the smaller barrier
  if E_reac > E_prod:
    E0 = E_reac
    dV1 = E_TS - E_reac
    dV2 = E_TS - E_prod
  else:
    E0 = E_prod
    dV1 = E_TS - E_prod
    dV2 = E_TS - E_reac

  # Convert units
  E0  *= cm_to_J  
  dV1 *= cm_to_J  
  dV2 *= cm_to_J
  vTS *= cm_to_J

  # Calculate intermediate terms to be placed into the kappa_e function 
  alpha1 = 2.0 * np.pi() * ( dV1 / (h*vTS) )
  alpha2 = 2.0 * np.pi() * ( dV2 / (h*vTS) )
  xi = ( (E-E0) / dV1 )

  a = ( 2.0*np.sqrt(alpha1*xi) ) / ( alpha1**(-1.0/2.0) + alpha2**(-1.0/2.0) )
  b = ( 2.0*np.sqrt( (xi-1.0)*alpha1+alpha2 ) ) / ( alpha1**(-1.0/2.0) + alpha2**(-1.0/2.0) )
  d = 2.0 * np.sqrt( alpha1*alpha2 - ( (4.0*np.pi()**2) / 16.0 ) )
  
  A = 2.0 * np.pi() * a
  B = 2.0 * np.pi() * b
  D = 2.0 * np.pi() * d

  # Calculate kappa(E)
  kappa_e = 1.0 - ( ( np.cosh(A-B) + np.cosh(D) ) / ( np.cosh(A+B) + np.cosh(D) ) )

  # Use a Fourier transform to change kappa(E) to kappa(T)
  kappa_t = kappa_e * sum()

  return kappa_t


def calc_wigner_tunn_coeff(vTS, T):
  ''' Calculates kappa(T) using the Wigner potential model. 
  '''
  
  kappa_T = 1.0 + (1.0/24.0) * ( (h*vTS) / (kB*T) )**2

  return kappa_T
