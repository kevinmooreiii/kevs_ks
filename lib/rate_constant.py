''' 
Calculates the rate constant
'''

import numpy as np
from constants import kB, h


def calc_canon_tst_rate_constant(Q_ts, Q_reac, barrier, kappa_T, T):
  ''' Calculates the Rate Constant using Canonical Transition State Theory. '''

  k_T = ( kappa_T * ((kB*T)/h) * (Q_ts/Q_reac) * np.exp(-barrier/ (kB*T) )

  return k_T

def calc_microcanon_tst_rate_constant(N_E, rho_E):
  ''' Calculates the Rate Constant using Microcanonical Transition State Theory. '''

  k_E = ( N_E / (h*rho_E) )  

  return k_E

