########################################################################
# Team Spectacular Stellars: Arian Andalib, Ashley Stone, Jonathan Kho, Emma Oswald
# AST 304, Fall 2022
# Michigan State University
########################################################################
"""
<This module is starter code to integrate stellar structure equations>
"""

import numpy as np
from eos import * # fill this in
from ode import rk4 # fill this in
from astro_const import * # fill this in

def oldpressure(rho, mue,K):
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
    
    p = K*(1/5)*(3/(8*(np.pi)))**(2/3)*(ac.h**2/ac.m_e)*(rho/(mue*ac.m_u))**(5/3)
    
    return p

def olddensity(p, mue):
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
    
    rho = ((p/((1/5)*(3/(8*(np.pi)))**(2/3)*(ac.h**2/ac.m_e)))**(3/5))*(mue*ac.m_u)
    
    return rho

def oldstellar_derivatives(m,z,mue):
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
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm
    """
    
    rho =olddensity(z[1], mue)
    
    dzdm = np.zeros_like(z)

    # evaluate dzdm
    drdm = (4*np.pi*z[0]**2*rho)**(-1)
    dPdm = (-G*m)/(4*np.pi*z[0]**4)
    
    dzdm = np.array([drdm, dPdm])
    
    return dzdm

def oldcentral_values(Pc,delta_m,mue):
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
    
    # compute initial values of z = [ r, p ]
    
    m = delta_m
    
    P = Pc
    rho = olddensity(P, mue)
    r = ((3*m)/(4*np.pi*rho))**(1/3)
    
    z = [r, P]
    
    return z
    
def oldlengthscales(m,z,mue):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg)
        z (array)
           [ r, p ] (units =[ m, Pascal])
        mue
            mean electron weight
    
    Returns
        z/|dzdm| (units =[ m, Pascal ])
    """

    # fill this in
    
    Hr = 4*np.pi*z[0]**3*olddensity(z[1],mue)
    
    Hp = 4*np.pi*z[0]**4*z[1]/(G*m)
    
    return np.array([Hr,Hp])
    
def oldintegrate(Pc,delta_m,eta,xi,mue,max_steps=10000):
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
    
    # set starting conditions using central values
    m = delta_m
    z = oldcentral_values(Pc,delta_m,mue)
    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        
        # are we at the surface?
        if (pressure < eta*Pc):
            break
            
        # store the step
        m_step[step] = m
        r_step[step] = z[0]
        p_step[step] = z[1]
        
        # set the stepsize
        h = xi*min(oldlengthscales(m,z,mue))
        
        # take a step
        z = rk4(oldstellar_derivatives,m,z,h,mue)
        m += h
        
        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]

def oldpressure_guess(m,mue,k):
    """
    Computes a guess for the central pressure based on virial theorem and mass-
    radius relation. 
    
    Arguments
        m
            mass of white dwarf (units are kg)
        mue
            mean electron mass
        
        k
            polytropic constant
    
    Returns
        P
            guess for pressure
    """
    # fill this in
    
    # For optional calibration
    # alpha = 
    # beta = 
    # Pguess = (alpha**5/beta**(20/3))*(G**5/Ke**4)*(m*mue**2)**(10/3)
    
    Pguess = (1/40.09269)*(G**5/k**4)*(m*mue**2)**(10/3)
    
    return Pguess


def Desired_mass(P_c, m, delta_m, eta, xi, mue):
    '''
    Finds the difference between the final mass found 
    by the integration and the specified desired mass
    
    Arguments
         P_c
            central pressure (units = Pascal)
         m
             desired mass
        delta_m
            initial offset from center (units = kg)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
    Returns
        m - m_step[-1]
    
    '''
    
    m_step, r_step, p_step = oldintegrate(P_c,delta_m,eta,xi,mue)

    return m - m_step[-1]
