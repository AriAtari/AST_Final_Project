""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

Team Jimmy Neutron: Anita Agasaveeran, Nick Persha, Levi Webb, Evelyn Bruinsma
"""

import numpy as np
from eos_template import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi
from reactions_template import pp_rate
from zams_template import Teff, surface_luminosity


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
    # fill this in
    Pc = 0.77 * G * (Msun*m)**2 / (Rsun*r)**4
    rhoc = 5.99 * (3 * m*Msun) / (fourpi * (Rsun*r)**3)
    Tc = 0.54 * (mu * m_u / kB) * (G * Msun*m / (Rsun*r))
    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
def stellar_derivatives(m,z, XH, pp_factor, central):
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure)
        XH
            mass fraction of ionized plasma
        pp_factor
            multiplicative factor for rate
        central (list)
            values of (P_c, rho_c, T_c)
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm, dL/dm
    """
    
    
    r = z[0]
    p = z[1]
    L = z[2]
    dzdm = np.zeros_like(z)
    
    rho, T = get_rho_and_T(p, central[0], central[1], central[2])

    # evaluate dzdm
    
    dzdm[0] = (4 * np.pi * r**2 * rho)**(-1)
    dzdm[1] = (-G * m)/(4 * np.pi * r**4)
    dzdm[2] = pp_rate(T,rho,pp_factor,XH)
    
    return dzdm

def central_values(central,delta_m,XH,pp_factor):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
        Pc
            central pressure (units = Pa) 
        delta_m
            core mass (units = kg) 
        XH
            mass fraction of ionized plasma
        pp_factor
            multiplicative factor for rate
        central (list)
            values of (P_c, rho_c, T_c)
    
    Returns
        z = array([ r, p, L ])
            central values of radius and pressure (units = kg, Pa, W) 
    """
    z = np.zeros(3)
    
    # compute initial values of z = [ r, p ] 
    
    r = ( (3 * delta_m) / (4 * np.pi * central[1]) )**(1/3)
    p = central[0]
    L = pp_rate(central[2],central[1],XH,pp_factor)*delta_m
    
    z[0] = r
    z[1] = p
    z[2] = L
    
    return z

def lengthscales(m,z,XH, pp_factor, central):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg) 
        z (array)
           [ r, p , L] (units = m, Pa, W) 
        XH
            mass fraction of ionized plasma
        pp_factor
            multiplicative factor for rate
        central (list)
            values of (P_c, rho_c, T_c)
    
    Returns
        z/|dzdm| (units = kg, kg)
    """

    # fill this in
    r = z[0]
    p = z[1]
    L = z[2]
    rho_hr, T = get_rho_and_T(p, central[0], central[1], central[2]) 
    
    Hr = 4 * np.pi * r**3 * rho_hr
    Hp = (4 * np.pi * r**4 * p) / (G * m)
    HL = L / pp_rate(T,rho_hr,XH,pp_factor)
   
    zdzdm = np.array([Hr,Hp])
    
    return zdzdm

def integrate(delta_m,eta,xi,central, XH, max_steps=10000, pp_factor = 1, args = ()):
    """
    Integrates the scaled stellar structure equations

    Arguments
        Pc
            central pressure (units = Pa)
        delta_m
            initial offset from center (units = kg)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        XH
            mass fraction of ionized plasma
        pp_factor
            multiplicative factor for rate
        central (list)
            values of (P_c, rho_c, T_c)
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (kg, m, Pa)
    """
        
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    L_step = np.zeros(max_steps)
    
    # set starting conditions using central values
    z = central_values(central,delta_m,XH,pp_factor)
    m = delta_m

    # stepsizes = lengthscales(delta_m,z,mue)
    # h = np.min(stepsizes) * xi

    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        luminosity = z[2]
        # are we at the surface?
        if (pressure < eta*central[0]):
            break
            
     
        # store the step       
        m_step[step] = m
        r_step[step] = radius
        p_step[step] = pressure
        L_step[step] = luminosity

        # set the stepsize
        stepsizes = lengthscales(m,z,XH, pp_factor, central)
        h = np.min(stepsizes) * xi
        
        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=(XH, pp_factor, central))
        m += h
        
        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps], L_step[0:Nsteps]

def Ldiff(Rfinal, Mwant, delta_m, eta, xi, XH, mu, pp_factor):
    thermal_c = central_thermal(Mwant, Rfinal, mu)
    L_nuc = integrate(delta_m, eta, xi, thermal_c, XH, args = (pp_factor))[3][-1]

    Rfinal = Rfinal * Rsun
    T_eff = Teff(Mwant)
    L_surface = surface_luminosity(T_eff, Rfinal)

    return L_nuc - L_surface
