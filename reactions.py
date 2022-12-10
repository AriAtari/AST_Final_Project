# 
#Team Spectacular Stellars: Arian Andalib, Ashley Stone, Jonathan Kho, Emma Oswald
# 

import numpy as np

def pp_rate(T,rho,XH= 0.7,pp_factor=1.0):
    """
    Specific heating rate from pp chain hydrogen burning. Approximate rate 
    taken from Hansen, Kawaler, & Trimble.
    
    Arguments
        T, rho
            temperature [K] and density [kg/m**3]
        XH
            mass fraction of hydrogen
        pp_factor
            multiplicative factor for rate
    Returns
        heating rate from the pp-reaction chain [W/kg]
    """
    
    # fill this out
    rate = pp_factor*((2.4e-3*rho*XH**2)/(T*1e-9)**(2/3)*np.exp(-3.380/(T*1e-9)**(1/3)))
    return rate

