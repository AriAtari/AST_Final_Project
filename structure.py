""" 
Routines for computing structure of fully convective star of a given mass and 
radius.
Team Spectacular Stellars: Arian Andalib, Ashley Stone, Jonathan Kho, Emma Oswald
"""

import numpy as np
from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi,Ke, sigmaSB
from reactions import pp_rate
from zams import Teff, surface_luminosity

def central_thermal(m,r,mu):
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
    R = r*Rsun
    Pc = 0.77*G*M**2/R**4
    rhoc = 5.99*(3*M)/(fourpi*R**3)
    Tc = 0.54*((mu*m_u)/kB)*(G*M)/R
    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
def stellar_derivatives(m,z,P_c,rho_c,T_c,pp_factor): # mu?
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure, luminsity)
        
        P_c, rho_c, T_c
            central pressure, density, and temperature in solar units
        
        pp_factor 
            multiplicative factor for rate
            
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm, dL/dm
    """
    rho, T = get_rho_and_T(z[1],P_c,rho_c,T_c)
        
    dzdm = np.zeros_like(z)

    # evaluate dzdm
    drdm = (4*np.pi*z[0]**2*rho)**(-1)
    dPdm = (-G*m)/(4*np.pi*z[0]**4)
    dLdm = pp_rate(T, rho,pp_factor = pp_factor)
    dzdm = np.array([drdm, dPdm, dLdm])
    
    return dzdm

def central_values(M,R,delta_m,mu, pp_factor=1.0):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
    
        M
            solar mass fraction
        R
            solar radius fraction
        delta_m
            core mass (units = kg)
        mu
            mean molecular weight
        pp_factor 
            multiplicative factor for rate
    
    Returns
        M
            solar mass fraction
        z
            initial values of (radius, pressure, luminsity)
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """
     
    # compute initial values of z = [ r, p, L ]
    Pc,rhoc,Tc = central_thermal(M,R,mu) # total mass from table
    
    r = ((3*Msun*delta_m)/(4*np.pi*rhoc))**(1/3)
    
    L = Msun*delta_m*pp_rate(Tc,rhoc,pp_factor = pp_factor)
    
    z = [r, Pc, L]
   
    return M, z, Pc, rhoc, Tc
    
def lengthscales(m,z,Pc,rhoc,Tc,pp_factor=1.0):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg)
            
        z (array)
           [ r, p, L ] (units =[ m, Pascal])
           
        P_c, rho_c, T_c
            central pressure, density, and temperature in solar units
            
        pp_factor 
            multiplicative factor for rate
        
    
    Returns
        z/|dzdm| = [Hr,Hp,HL] (units =[ m, Pascal ])
    """

    # fill this in
    
    rho, T = get_rho_and_T(z[1],Pc,rhoc,Tc)
        
    Hr = 4*np.pi*z[0]**3*rho
    
    Hp = 4*np.pi*z[0]**4*z[1]/(G*m)
    
    HL = 4*np.pi*z[0]**2*sigmaSB*T**4 
    
    return np.array([Hr,Hp,HL])
    
def integrate(M,r,delta_m,eta,xi,comp,max_steps=10000,pp_factor=1.0):
    """
    Integrates the scaled stellar structure equations
    Arguments
        M
            current mass coordinate (units = kg)
        r
            current radius coordinate (units = m)

        delta_m
            initial offset from center (units = kg)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        comp
            a 2D array of the chemical makeup of the elements in the star which 
            consists of 1D arrays for Z, A, and X values.    ###### constant throughout the star
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
        pp_factor 
            multiplicative factor for rate                       
               
    Returns
        m_step, r_step, p_step, L_step
            arrays containing mass coordinates, radii, pressures, and luminosities during 
            integration (units:[kg,m,Pascal,W])
    """
    
    z_atom = comp[0]
    a = comp[1]
    x = comp[2]
    
    mu = mean_molecular_weight(z_atom,a,x) ##### calc once
    
    M, z, Pc, rhoc, Tc = central_values(M,r,delta_m, mu, pp_factor)
     
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
        if (z[1] < eta*Pc):
            break
            
        # store the step
        m_step[step] = m
        r_step[step] = z[0]
        p_step[step] = z[1]
        L_step[step] = z[2]
        
        # set the stepsize
        h = xi*min(lengthscales(m,z,z[1],rhoc,Tc,pp_factor)) ################ mass coordinate (little m)
        
        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=(Pc,rhoc,Tc,pp_factor)) #################### little m,
        m += h
        
        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps],L_step[0:Nsteps]


