""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

<Team name, members go here>
"""

import numpy as np
from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi
from reactions import pp_rate

def central_thermal(m,r,mu):
    """ 
    Computes the central pressure, density, and temperature from the polytropic
    relations for n = 3/2.

    Arguments
        m
            mass in solar units
        r
            radius is solar units
        mu
            mean molecular weight
    Returns
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """

    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
def stellar_derivatives(m,z,mue,rate):
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure)
        mue
            ratio, nucleons to electrons.  For a carbon-oxygen white dwarf, 
            mue = 2.
        rate
            the nuclear heating rate per mass unit calculated using
            the reaction routine pp_rate.
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm, dL/dm
    """
    
    rho = density(z[1], mue)
    
    dzdm = np.zeros_like(z)

    # evaluate dzdm
    drdm = (4*np.pi*z[0]**2*rho)**(-1)
    dPdm = (-G*m)/(4*np.pi*z[0]**4)
    dLdm = rate
    
    dzdm = np.array([drdm, dPdm, dLdm])
    
    return dzdm

def central_values(Pc,delta_m,mue,mu):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
        Pc
            central pressure (units = Pascal)
        delta_m
            core mass (units = kg)
        mue
            nucleon/electron ratio
    
    Returns
        z = array([ r, p ])
            central values of radius and pressure (units =[ m, Pascal])
    """
    
    # compute initial values of z = [ r, p, L ]
    
    m = delta_m
    
    P = Pc
    rho = density(P, mue)
    r = ((3*m)/(4*np.pi*rho))**(1/3)
    
    central = central_therm(m,r,mu)
    Tc = central[2]
    rhoc = central[1]
    
    L = delta_m*pp_rate(Tc,rhoc)
    
    z = [r, P, L]
    
    return z
    
def lengthscales(m,z,mue,R,Teff):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg)
        z (array)
           [ r, p, L ] (units =[ m, Pascal])
        mue
            mean electron weight
        R 
            outer radius of star
        Teff
            the calculated effective temperature
    
    Returns
        z/|dzdm| (units =[ m, Pascal ])
    """

    # fill this in
    
    Hr = 4*np.pi*z[0]**3*density(z[1],mue)
    
    Hp = 4*np.pi*z[0]**4*z[1]/(G*m)
    
    HL = 4*np.pi*R**2*sigmaSB*Teff**4
    
    return np.array([Hr,Hp,HL])
    
def integrate(Pc,delta_m,eta,xi,mue,max_steps=10000):
    """
    Integrates the scaled stellar structure equations

    Arguments
        Pc
            central pressure (units = Pascal)
        delta_m
            initial offset from center (units = kg)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (units:[kg,m, Pascal])
    """
        
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    L_step = np.zeros(max_steps)
    
    # set starting conditions using central values
    m = delta_m
    z = central_values(Pc,delta_m,mue)
    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        luminosity = z[2]
        
        # are we at the surface?
        if (pressure < eta*Pc):
            break
            
        # store the step
        m_step[step] = m
        r_step[step] = z[0]
        p_step[step] = z[1]
        L_step[step] = z[2]
        
        # set the stepsize
        h = xi*min(lengthscales(m,z,mue,R,Teff))
        
        # take a step
        z = rk4(stellar_derivatives,m,z,h,mue,R,Teff)
        m += h
        
        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps],L_step[0:Nsteps]