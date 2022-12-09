""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

<Team name, members go here>
"""

import numpy as np
<<<<<<< HEAD
from eos import *
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi, Ke
from reactions import pp_rate

def central_thermal(m,r,mu):
=======
<<<<<<<< HEAD:structure.py
from eos import get_rho_and_T, mean_molecular_weight, density, pressure
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi,Ke
========
from eos import *
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi, Ke
>>>>>>>> origin/Emma:structure_template.py
from reactions import pp_rate
from zams import Teff, surface_luminosity

def central_thermal(m,R,mu):
>>>>>>> origin/Emma
    """ 
    Computes the central pressure, density, and temperature from the polytropic
    relations for n = 3/2.

    Arguments
        m
<<<<<<< HEAD
            mass in solar units
=======
            mass in solar units  total mass
>>>>>>> origin/Emma
        r
            radius is solar units
        mu
            mean molecular weight
    Returns
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """
    # fill this in
<<<<<<< HEAD
    Pc = 0.0
    rhoc = 0.0
    Tc = 0.0
=======
    M=m*Msun
    Pc = 0.77*G*M**2/R**4
    rhoc = 5.99*(3*M)/(fourpi*R**3)
    Tc = 0.54*((mu*m_u)/kB)*(G*M)/R
>>>>>>> origin/Emma
    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
<<<<<<< HEAD

def stellar_derivatives(m,z,mue):
=======
<<<<<<<< HEAD:structure.py
def stellar_derivatives(m,z,Teff,rho,mue=2):
========

def stellar_derivatives(m,z,mue):
>>>>>>>> origin/Emma:structure_template.py
>>>>>>> origin/Emma
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
<<<<<<< HEAD
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm
=======
        rate
            the nuclear heating rate per mass unit calculated using
            the reaction routine pp_rate.
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm, dL/dm
>>>>>>> origin/Emma
    """
    
    rho = density(z[1], mue)
    
    dzdm = np.zeros_like(z)

    # evaluate dzdm
    drdm = (4*np.pi*z[0]**2*rho)**(-1)
    dPdm = (-G*m)/(4*np.pi*z[0]**4)
<<<<<<< HEAD
    
    dzdm = np.array([drdm, dPdm])
    
    return dzdm

def central_values(m,r,delta_m,pp_factor):
=======
    dLdm = pp_rate(Teff, rho)
    
    dzdm = np.array([drdm, dPdm, dLdm])
    
    return dzdm

def central_values(M,r,delta_m, mu, pp_factor=1.0, mue=2): ##############################################
>>>>>>> origin/Emma
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
<<<<<<< HEAD
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
    rho = density(P, mue)
    r = ((3*m)/(4*np.pi*rho))**(1/3)
    
    z = [r, P]
    
    return z
    
def lengthscales(m,z,mue):
=======

        mue
            nucleon/electron ratio ?????????
            
        M
            solar mass fraction
        r
            total radius given by table???????
        delta_m
            core mass (units = kg)
        pp_factor 
            
    
    Returns
        z = array([ r, pc ])
            central values of radius and pressure (units =[ m, Pascal])
    """
     
    # compute initial values of z = [ r, p, L ]
    central = central_thermal(M,r,mu) ####################################### total mass from table
    Tc = central[2]
    rhoc = central[1]
    Pc = central[0]

    P = Pc
    
    rho = density(Pc, mue) 
    
    r = ((3*delta_m)/(4*np.pi*rho))**(1/3)
    
    L = delta_m*pp_rate(Tc,rhoc,pp_factor)
    
    z = [r, P, L]
   
    return M, z, Pc, rhoc, Tc
    
def lengthscales(m,z,Pc,rhoc,T_eff,pp_factor=1.0):
>>>>>>> origin/Emma
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg)
        z (array)
<<<<<<< HEAD
           [ r, p ] (units =[ m, Pascal])
        mue
            mean electron weight
=======
           [ r, p, L ] (units =[ m, Pascal])
        R 
            outer radius of star
        Teff
            the calculated effective temperature
>>>>>>> origin/Emma
    
    Returns
        z/|dzdm| (units =[ m, Pascal ])
    """

    # fill this in
    
    Hr = 4*np.pi*z[0]**3*density(z[1],mue)
    
    Hp = 4*np.pi*z[0]**4*z[1]/(G*m)
    
<<<<<<< HEAD
    return np.array([Hr,Hp])
    
def integrate(m,r,delta_m,eta,xi,pp_factor,max_steps=10000):
=======
    HL = 4*np.pi*z[0]**2*sigmaSB*T_eff**4 ##############################################
    
    return np.array([Hr,Hp,HL])
    
def integrate(M,r,delta_m,eta,xi,comp,max_steps=10000,pp_factor=1.0):
>>>>>>> origin/Emma
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
<<<<<<< HEAD
=======
        comp
            a 2D array of the chemical makeup of the elements in the star which 
            consists of 1D arrays for Z, A, and X values.    ###### constant throughout the star
>>>>>>> origin/Emma
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (units:[kg,m, Pascal])
    """
<<<<<<< HEAD
        
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    
    # set starting conditions using central values
    m = delta_m
    z = central_values(Pc,delta_m,mue)
    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
=======
    
    z_atom = comp[0]
    a = comp[1]
    x = comp[2]
    
    mu = mean_molecular_weight(z_atom,a,x) ##### calc once
    
    M, z, rhoc, Tc = central_values(M,r,delta_m, mu, pp_factor=1.0)
     
    T_eff = Teff(M)
    
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    L_step = np.zeros(max_steps)
    
    # set starting conditions using central values
    m = delta_m
    
    Nsteps = 0
    for step in range(max_steps):
#         radius = central_vals[0]
#         pressure = central_vals[1]
#         luminosity = central_vals[2]           don't need?
#         rho = density(pressure, mue)
        
>>>>>>> origin/Emma
        
        # are we at the surface?
        if (pressure < eta*Pc):
            break
            
        # store the step
        m_step[step] = m
        r_step[step] = z[0]
        p_step[step] = z[1]
<<<<<<< HEAD
        
        # set the stepsize
        h = xi*min(lengthscales(m,z,mue))
        
        # take a step
        z = rk4(stellar_derivatives,m,z,h,mue)
=======
        L_step[step] = z[2]
        
        # set the stepsize
        h = xi*min(lengthscales(m,z,p[step],rhoc,T_eff,pp_factor=1.0)) ########################################### mass coordinate (little m)
        
        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=[T_eff,rhoc]) ################################################## little m
>>>>>>> origin/Emma
        m += h
        
        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
<<<<<<< HEAD
=======
<<<<<<<< HEAD:structure.py
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps],L_step[0:Nsteps]









========
>>>>>>> origin/Emma
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]

def pressure_guess(m,mue):
    """
    Computes a guess for the central pressure based on virial theorem and mass-
    radius relation. 
    
    Arguments
        m
            mass of white dwarf (units are kg)
        mue
            mean electron mass
    
    Returns
        P
            guess for pressure
    """
    # fill this in
    
    # For optional calibration
    # alpha = 
    # beta = 
    # Pguess = (alpha**5/beta**(20/3))*(G**5/Ke**4)*(m*mue**2)**(10/3)
    
<<<<<<< HEAD
    Pguess = (1/40.09269)*(G**5/Ke**4)*(m*mue**2)**(10/3)
    
    return Pguess
=======
    Pguess = (1/40.09269)*(G**5/Ke**4)*(m*mue**2)**(10/3) # (1/40.09269)*
    
    return Pguess
>>>>>>>> origin/Emma:structure_template.py
>>>>>>> origin/Emma
