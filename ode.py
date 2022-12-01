########################################################################
# Team Spectacular Stellars: Arian Andalib, Ashley Stone, Jonathan Kho, Emma Oswald
# AST 304, Fall 2020
# Michigan State University
########################################################################


"""
This module contains functions for approximating a positional value for a given ordinary differential equation. The functions are
    - fEuler (Forward Euler Approximation)
    - rk2 (Second-order Runge-Kutta)
    - rk4 (Fourth-order Runge-Kutta)
"""

# all routines that take a single step should have the same interface
# fEuler is complete, except for documentation. you can use this as a pattern 
# for the other two routines.
def fEuler(f,t,z,h,args=()):
    """
    Routine for forward Euler approximation. A first-order approximation to solve a given
    ODE. Called with one function, a starting time value, a starting z value, a time step, 
    and any additional arguments to pass to the inputed function. Returns an estimated z 
    (position) value for the given ODE and initial conditions.
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        t (float)
            the initial time value
            
        z (float)
            initial z position
            
        h (float)
            time step value
    
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
    
    # The following trick allows us to pass additional parameters to f
    # first we make sure that args is of type tuple; if not, we make it into
    # that form
    if not isinstance(args,tuple):
        args = (args,)
    
    # when we call f, we use *args to pass it as a list of parameters.
    # for example, if elsewhere we define f like
    # def f(t,z,x,y):
    #    ...
    # then we would call this routine as
    # znew = fEuler(f,t,z,h,args=(x,y))
    #
    return z + h*f(t,z,*args)

# You will need to flesh out the following routines for a second-order
# Runge-Kutta step and a fourth order Runge-Kutta step.

def rk2(f,t,z,h,args=()):
    """
    Routine for a second order Runge-Kutta step. A second-order approximation to solve a
    givenODE. Called with one function, a starting time value, a starting z value, a time 
    step, and any additional arguments to pass to the inputed function. Returns an  
    estimated z (position) value for the given ODE and initial conditions.
    
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        t (float)
            the initial time value
            
        z (float)
            initial z position
            
        h (float)
            time step value
    
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
    
    if not isinstance(args,tuple):
        args = (args,)
    
    #estimate of the midpoint of the interval
    z_p = (z+(h/2)*f(t, z,*args))

    #calculating z(t+h) using z_p
    z_new = (z + h*f((t+(h/2)), z_p,*args))
    
    return z_new

def rk4(f,t,z,h,args=()):
    """
    Routine for a fourth order Runge-Kutta step. A fourth-order approximation to solve a
    givenODE. Called with one function, a starting time value, a starting z value, a time 
    step, and any additional arguments to pass to the inputed function. Returns an 
    estimated z (position) value for the given ODE and initial conditions.
    
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        t (float)
            the initial time value
            
        z (float)
            initial z position
            
        h (float)
            time step value
            
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """ 
   
    if not isinstance(args,tuple):
        args = (args,)
    
    #four estimates for dz/dt
    k1 = f(t, z,*args)
    k2 = f(t+(h/2), z+(h/2)*k1,*args)
    k3 = f(t+(h/2), z+(h/2)*k2,*args)
    k4 = f(t+h, z+h*k3,*args)
    
    #calculate the weighted avg. to find z(t+h)
    znew = z +(h/6)*(k1+2*k2+2*k3+k4)
 
    return znew