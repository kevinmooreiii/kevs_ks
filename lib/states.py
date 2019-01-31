"""
  Calculates information about the number and density of states used in Microcanonical TST
  Supported: ONLY Vibrational States
"""

from constants import cm_to_J


def calc_sum_of_states(freqs, E0, E_max):
    """ Use the Beyer-Swinehart Algorithm to compute Sum-of-States (N_E) up to the requested energy level (E).
        N_E is a list where the index of the list indicates the energy level.
        For example: N_E[1435] = 10^12 implies there are 10^12 vibrational states at 1435 *units* level
    """

    # Convert E0 from to cm-1
    E0 = E0 / cm_to_J

    # Initialize a list for for energies greater than TS energy, E0
    # All states below E0 set to zero since they do not contribute
    N_0 = [0 for E in range(0, E0+1)]
    N_1 = [1 for E in range(E0+1, E_max+1)]
    N_E = N_0 + N_1

    # Loop over the freqs and count up the number of freqs
    for v in freqs:
        for n in range(E0+1+v, Emax+1):
            N_E[n] = N_E[n] + N_E[n-v]

    return N_E


def calc_density_of_states(freqs, E_max):
    """ Determine the density of vibrational states: pE, which is initialized by assigning 1 state per 1 cm-1 """

    # Initialize list for the density states
    p_E = [0 for E in range(E_max+1)]
    p_E[0] = 1

    # Compute rho_E by finite differences of the sum-of-states
    for v in freqs:
        for n in range(v, Emax+1):
            p_E[n] = p_E[n] + p_E[n-v]

    return p_E

