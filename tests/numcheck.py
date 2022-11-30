import numpy as np

def within_tolerance(v,a,tol=None):
    if not tol:
        tol = 4.0*np.finfo(1.0).eps
    
    return np.abs(1.0-v/a) < tol
