'''
  Calculates several energy quantities
'''

from constants import cm_to_J, hart_to_J


def calc_H0(E_elec, freqs):
  ''' Calculates the Enthalpy at 0 K. '''

  # Calculate the zero-point vibrational energy (ZPVE) wihin the Harmonic-Oscillator approximation
  E_zpve = calc_zpve(freqs) * cm_to_J

  # Add the ZPVE to the electronic energy to get the Enthalpy at 0 K (H_0), while converting to Joules
  H_0 = ( E_elec * hart_to_J) + ( E_zpve * cm_to_J )

  return H_0


def calc_zpve(freqs):
  ''' Calculates the zero-point vibrational energy. '''

  E_zpve = (1.0/2.0) * sum(freqs)

  return E_zpve

