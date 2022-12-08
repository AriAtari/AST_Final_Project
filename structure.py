""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

<Team name, members go here>
"""

import numpy as np
from eos import get_rho_and_T, mean_molecular_weight, density, pressure
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi,Ke
from reactions import pp_rate
from zams import Teff, surface_luminosity

def central_thermal(m,R,mu):
    """ 
    Computes the central pressure, density, and temperature from the polytropic
    relations for n = 3/2.

    Arguments
        m
            mass in solar units  total mass
        r
            radius is solar units
        mu
            mean molecular weight
    Returns
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """
    # fill this in
    M=m*Msun
    Pc = 0.77*G*M**2/R**4
    rhoc = 5.99*(3*M)/(fourpi*R**3)
    Tc = 0.54*((mu*m_u)/kB)*(G*M)/R
    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
def stellar_derivatives(m,z,Teff,rho,mue=2):
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
    dLdm = pp_rate(Teff, rho)
    
    dzdm = np.array([drdm, dPdm, dLdm])
    
    return dzdm

def central_values(M,r,delta_m, mu, pp_factor=1.0, mue=2): ##############################################
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments

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
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg)
        z (array)
           [ r, p, L ] (units =[ m, Pascal])
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
    
    HL = 4*np.pi*z[0]**2*sigmaSB*T_eff**4 ##############################################
    
    return np.array([Hr,Hp,HL])
    
def integrate(M,r,delta_m,eta,xi,comp,max_steps=10000,pp_factor=1.0):
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
        comp
            a 2D array of the chemical makeup of the elements in the star which 
            consists of 1D arrays for Z, A, and X values.    ###### constant throughout the star
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (units:[kg,m, Pascal])
    """
    
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
        
        
        # are we at the surface?
        if (pressure < eta*Pc):
            break
            
        # store the step
        m_step[step] = m
        r_step[step] = z[0]
        p_step[step] = z[1]
        L_step[step] = z[2]
        
        # set the stepsize
        h = xi*min(lengthscales(m,z,p[step],rhoc,T_eff,pp_factor=1.0)) ########################################### mass coordinate (little m)
        
        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=[T_eff,rhoc]) ################################################## little m
        m += h
        
        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps],L_step[0:Nsteps]









