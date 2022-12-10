
Contents
--------

0. `README.md`: this file
1. `astro_const.py`: module containing physical constants, uses `astropy`
2. `eos.py`: Routines to compute an adiabatic equation of state.
3. `structure.py`: starter code to integrate stellar structure equations
4. `ode.py` : This module contains functions for approximating a positional value for a given ordinary differential equation. The functions are
    - fEuler (Forward Euler Approximation)
    - rk2 (Second-order Runge-Kutta)
    - rk4 (Fourth-order Runge-Kutta)
5. `reactions.py` : contatins pp_rate() which is Specific heating rate from pp chain hydrogen burning. Approximate rate 
    taken from Hansen, Kawaler, & Trimble.
6. `zamz.py` : Routines for computing the zero-aged main-sequence for low-mass stars.
7. `old_eos.py` : Routines from computational project 2 to complete question 1 (confirming scalings).

More Documentation!
--------

- `eos.py`:

> mean_molecular_weight(Z,A,X):
	Computes the mean molecular weight for a fully ionized plasma with an arbitrary mixture of species

    Arguments
        Z, A, X (either scaler or array)
            charge numbers, atomic numbers, and mass fractions
            The mass fractions must sum to 1

	Returns
		mu
   
> get_rho_and_T(P,P_c,rho_c,T_c):
    
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


--------------------------------------------------
- `structure.py`:

> central_thermal(m,r,mu):
    
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

 
> stellar_derivatives(m,z,P_c,rho_c,T_c): 
    
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure, luminsity)
        
	  P_c, rho_c, T_c
            central pressure, density, and temperature in solar units
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm, dL/dm
 
> central_values(M,R,delta_m,mu, pp_factor=1.0): 
    
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
    
        M
            solar mass fraction
        R
            solar radius fraction
        delta_m
            change in M (units = kg)
        pp_factor 
            multiplicative factor for rate
    
    Returns
        z = array([ r, pc ])
            central values of radius and pressure (units =[ m, Pascal])
    
          
> lengthscales(m,z,P_c,rho_c,T_c,pp_factor=1.0):
    
   Computes the radial length scale H_r and the pressure length H_P and luminosity length scale H_L.
    
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
        z/|dzdm| (units =[ m, Pascal, Watts ])

> integrate(M,r,delta_m,eta,xi,comp,max_steps=10000,pp_factor=1.0):
    
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
        mue
            mean electron mass
        comp
            a 2D array of the chemical makeup of the elements in the star which 
            consists of 1D arrays for Z, A, and X values.    ###### constant throughout the star
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
	  pp_factor 
            multiplicative factor for rate 	                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (units:[kg,m, Pascal])


-------------------------------------------

- `zams.py`

> Teff(Mwant):
    
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
    
> surface_luminosity(Teff,radius):
    
    Photospheric luminosity [W]
    
    Arguments
        Teff [K]
        radius [m]
    