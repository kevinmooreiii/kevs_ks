"""
  Calculates several thermodynamic energies
"""

from constants import cm_to_J, hart_to_J


def calc_enthalpy_0K(E_elec, freqs):
    """ Calculates the Enthalpy at 0 K. """

    # Calculate the zero-point vibrational energy (ZPVE) wihin the Harmonic-Oscillator approximation
    E_zpve = calc_zpve(freqs)

    # Add the ZPVE to the electronic energy to get the Enthalpy at 0 K (H_0), while converting to Joules
    H_0K = (E_elec * hart_to_J) + (E_zpve * cm_to_J)

    return H_0K


def calc_zpve(freqs):
    """ Calculates the zero-point vibrational energy. """

    E_zpve = (1.0/2.0) * sum(freqs)

    return E_zpve


def calc_enthalpy_TK(H_0K):
    """ Calculates the total enthalpy using the Rigid-Rotor Harmonic Oscillator approximations.
    """

    # Calc contributions from vibrational, translational, rotation, and electronic motion
    H_vib = calc_vib_contrib_enthalpy()
    H_trans = calc_vib_contrib_enthalpy()
    H_rot = calc_vib_contrib_enthalpy()
    H_elec = calc_vib_contrib_enthalpy()

    # Calculate enthalpy at temperature T
    H_TK = H_0K + H_vib + H_trans + H_rot + H_elec

    return H_TK


def calc_entropy_TK(H_0K):
    """ Calculates the total entropy using the Rigid-Rotor Harmonic Oscillator approximations.
    """

    # Calc contributions from vibrational, translational, rotation, and electronic motion
    S_vib = calc_vib_contrib_entropy()
    S_trans = calc_vib_contrib_entropy()
    S_rot = calc_vib_contrib_entropy()
    S_elec = calc_vib_contrib_entropy()

    # Calculate enthalpy at temperature T
    S_TK = S_vib + S_trans + S_rot + S_elec

    return S_TK


def calc_Gibbs_TK(T, H_TK, S_TK)
    """ Calculates the Gibb's Free Energy (G_TK) at temperature T using the enthalpy and entropy.
    """

    G_TK = H_TK - T*S_TK

    return G_TK
