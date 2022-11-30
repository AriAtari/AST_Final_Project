import pytest
from reactions import pp_rate
from numcheck import within_tolerance

def analytical_Texponent(T):
    """Analytical value of the logarithmic derivative of the pp heating rate
    (eq. [4] of the instructions)
    
    Arguments
        temperature [K]
    """
    return (-2 + 3.38/(T*1.0e-9)**(1/3))/3

def Texponent(T,XH,rho,pp_factor,delta=0.0001,):
    """Numerical estimate of the logarithmic derivative of the pp heating rate
    computed via central differencing.
    
    Arguments
        temperature [K]
        hydrogen mass fraction
        mass density [kg m**-3] 
        scaling factor
        fractional difference in temperature
    """
    Thigh = T*(1.0+delta)
    Tlow = T*(1.0-delta)
    r = pp_rate(T,rho,XH,pp_factor)
    rhigh = pp_rate(Thigh,rho,XH,pp_factor)
    rlow = pp_rate(Tlow,rho,XH,pp_factor)
    return T*(rhigh-rlow)/(Thigh-Tlow)/r

def test_pp():
    XH = 0.7
    rho = 1.0e3
    T = 1.0e7
    f = 1.0
    ref_value = 3.8930e-06
    base_rate = pp_rate(T,rho,XH,f)
    rate_XH = pp_rate(T,rho,XH/2,f)
    rate_rho = pp_rate(T,rho*10,XH,f)
    rate_f = pp_rate(T,rho,XH,f*10)
    
    # check against reference value
    assert within_tolerance(base_rate,ref_value,tol=1.0e-4),\
        "rate has incorrect value at T = 1.0e7 K, rho = 1000 kg m**-3, XH = 0.7"
    
    # check scalings with density, hydrogen mass fraction, scale factor
    assert within_tolerance(rate_XH,base_rate/4),\
        "rate does not scale properly with X_H"
    assert within_tolerance(rate_rho,base_rate*10),\
        "rate does not scale properly with rho"
    assert within_tolerance(rate_f,base_rate*10),\
        "rate does not scale properly with pp_factor"

    # check temperature scaling at T = 1.0e7 K and 2.0e7 K
    n_num = Texponent(T,XH,rho,f)
    n = analytical_Texponent(T)
    assert within_tolerance(n_num,n,tol=1.0e-3),\
        "rate has wrong temperature exponent"
    n_num = Texponent(2*T,XH,rho,f)
    n = analytical_Texponent(2*T)
    assert within_tolerance(n_num,n,tol=1.0e-3),\
        "rate has wrong temperature exponent"
    