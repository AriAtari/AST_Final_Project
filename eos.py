"""
Routines to compute an adiabatic equation of state.

<Team name, members go here>
"""

import astro_const as ac
import numpy as np
import astro_const as ac

def mean_molecular_weight(Z,A,X):
    """Computes the mean molecular weight for a fully ionized plasma with an
    arbitrary mixture of species

    Arguments
        Z, A, X (either scaler or array)
            charge numbers, atomic numbers, and mass fractions
            The mass fractions must sum to 1
    """
    Zs = np.array(Z)
    As = np.array(A)
    Xs = np.array(X)
    assert np.sum(Xs) == 1.0

    # compute value of mean molecular weight
    mu = np.sum(Xs*(1+Zs)/As)**(-1)
    return mu

def get_rho_and_T(P,P_c,rho_c,T_c):
    """
    Compute density and temperature along an adiabat of index gamma given a
    pressure and a reference point (central pressure, density, and temperature).

    Arguments
        P (either scalar or array-like)
            value of pressure

        P_c, rho_c, T_c
            reference values; these values should be consistent with an ideal
            gas EOS

    Returns
        density, temperature
    """

    # replace with computed values
    rho = rho_c*(P/P_c)**(3/5)
    T = T_c*(P/P_c)**(1-3/5)

    return rho, T

def pressure(rho, mue, K = 1): # Added polytropic constant K
    """
    Arguments
        rho
            mass density (kg/m**3)
        mue
            baryon/electron ratio

    Returns
        electron degeneracy pressure (Pascal)
    """

    # replace following lines with body of routine

    # p = K*rho**(5/3) #  Added polytropic constant K

    p = K*(1/5)*(3/(8*(np.pi)))**(2/3)*(ac.h**2/ac.m_e)*(rho/(mue*ac.m_u))**(5/3)

    return p

def density(p, mue, K = 1):
    """
    Arguments
        p
            electron degeneracy pressure (Pascal)
        mue
            baryon/electron ratio

    Returns
        mass density (kg/m**3)
    """

    # replace following lines with body of routine

    # rho = (p/K)**3/5

    rho = (K**(-3/5))*((p/((1/5)*(3/(8*(np.pi)))**(2/3)*(ac.h**2/ac.m_e)))**(3/5))*(mue*ac.m_u)

    return rho
