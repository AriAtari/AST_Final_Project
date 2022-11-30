import pytest
from astro_const import c,G,h,Msun,Rsun,Lsun
from numcheck import within_tolerance

def test_const():
    assert within_tolerance(c,2.99792458e+08)
    assert within_tolerance(G,6.67430e-11,tol=1.0e-6)
    assert within_tolerance(h,6.626e-34,tol=1.0e-4)
    assert within_tolerance(Msun,1.9884e+30,tol=1.0e-4)
    assert within_tolerance(Rsun,6.957e+08,tol=1.0e-4)
    assert within_tolerance(Lsun,3.828e+26,tol=1.0e-4)

    assert within_tolerance(c,2.99792458e+08)
    assert within_tolerance(G,6.67430e-11,tol=1.0e-6)
    assert within_tolerance(h,6.626e-34,tol=1.0e-4)
    assert within_tolerance(Msun,1.9884e+30,tol=1.0e-4)
    assert within_tolerance(Rsun,6.957e+08,tol=1.0e-4)
    assert within_tolerance(Lsun,3.828e+26,tol=1.0e-4)
