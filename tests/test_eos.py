import pytest
from eos import get_rho_and_T, mean_molecular_weight
from numcheck import within_tolerance

def test_chem():
    #scaling for H,He mixture
    for X in [0.0,0.3,0.7,1.0]:
        mu = mean_molecular_weight([1.0,2.0],[1.0,4.0],[X,1.0-X])
        assert within_tolerance(mu,4.0/(5.0*X+3.0))    

def test_adiabat():
    rhoc = 1.0e4
    Tc = 1.0e7
    eos_const=4.0
    Pc = eos_const*rhoc*Tc
    Ptest = Pc/10
    gamma = 5/3
    rho, T = get_rho_and_T(Ptest,Pc,rhoc,Tc)
    
    # is P/rho**gamma = const
    assert within_tolerance(Ptest/rho**gamma,Pc/rhoc**gamma),\
        "density does not follow adiabatic relation"
    # is EOS preserved?
    assert within_tolerance(Ptest/rho/T,eos_const),\
        "temperature does not follow adiabatic relation"
