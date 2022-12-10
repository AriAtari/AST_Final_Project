"""
Routines for computing the zero-aged main-sequence for low-mass stars.
Team Spectacular Stellars: Arian Andalib, Ashley Stone, Jonathan Kho, Emma Oswald
"""

import numpy as np
from astro_const import fourpi, sigmaSB, Rsun

def Teff(Mwant):
    """
    Interpolates effective temperatures given mass from tabulated [1] values 
    for low-mass stars.
    [1] Chabrier, Baraffe, Allard, and Hauschildt. Evolutionary Models for Very 
    Low-Mass Stars and Brown Dwarfs with Dusty Atmospheres. Astrophys. Jour. 
    542:464--472, Oct. 2000.
    Parameters
        Mwant (float, scalar or array)
            Mass of the star(s) in units of solar masses
    Returns
       Teff (float, same type/shape as Mwant)
            Interpolated effective temperatures of the stars in Kelvin.
    """

    # tabulated values from Chabrier et al. (2000)
    masses = np.array([0.1,0.15,0.2,0.3]) # [Msun]
    Teffs = np.array([2800.0,3150.0,3300.0,3400.0]) # [K]
    
    # fill this out to perform interpolation to find Teff for Mwant
    Teff = np.interp(Mwant,masses,Teffs)
    return Teff

def surface_luminosity(Teff,radius):
    """
    Photospheric luminosity [W]
    
    Arguments
        Teff [K]: Effective temperature of star
        radius [m]: Radius of star's surface
    Returns
        luminosity: Luminosity from star's surface
    """
    
    # fill this in
    luminosity = 4*np.pi*(radius)**2*sigmaSB*Teff**4
    return luminosity