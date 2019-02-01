"""
Libraries for computing rate constants
  SUPPORTED: Canoncial Transition State Theory, Microcanonical RRKM Theory
"""

import numpy as np
from constants import kB, h


def calc_canon_tst_rate_constant(T, Q_ts, Q_reac, barrier, kappa_T):
    """ Calculates the Rate Constant using 
        Canonical Transition State Theory. 
    """

    k_T =  kappa_T * ((kB*T)/h) * (Q_ts/Q_reac) * np.exp(-barrier/ (kB*T) 

    return k_T


def calc_rrkm_ke_rate_constant(N_E, p_E):
    """ Calculates the Rate Constant 
        using Microcanonical Transition State Theory. 
    """

    k_E =  N_E / (h*p_E)   

    return k_E
