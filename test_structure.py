import pytest
from structure import central_thermal
from numcheck import within_tolerance

def test_central_thermal():
    # check value for solar mass, radius and mean molecular weight 0.6
    Pcs,rhocs,Tcs = central_thermal(1.0,1.0,0.6)
    assert within_tolerance(Pcs,8.675e+14,tol=1.0e-2)
    assert within_tolerance(rhocs,8.445e+03,tol=1.0e-2)
    assert within_tolerance(Tcs,7.413e+06,tol=1.0e-2)
    
    # check scalings with mass, radius by trying 0.5 solar mass, 0.7 solar 
    # radius.
    Pc1, rhoc1, Tc1 = central_thermal(0.5,0.7,0.6)
    assert within_tolerance(Pc1, Pcs*0.5**2/0.7**4)
    assert within_tolerance(rhoc1,rhocs*0.5/0.7**3)
    assert within_tolerance(Tc1,Tcs*0.5/0.7)
        