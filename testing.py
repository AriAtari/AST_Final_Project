"""
This script runs the various testing functions for the ZAMS project. This script can be executed from the main project directory. It assumes that the various testing scripts exist in a subdirectory named `tests`. 

Usage:

When ready to test the various functions needed to compute the ZAMS structure of low-mass, fully-convective stars, simply uncomment the appropriate function calls below (see Section 2 of the Instructions). You do NOT need to delete the various lines from the testing functions that are referenced in the Instructions (this was out-of-date info...sorry). Then, from your main project directory containing the Python source files for your project, simply execute:

> python testing.py 

And the scripts should run as expected and return an error if your code fails the check. 

Note:

If you do not have the pytest package installed, you do NOT need it. Simply comment out the `import pytest` line in each of the testing source files. 

"""
import sys
sys.path.append('tests')
from test_eos import *
from test_reactions import *
from test_structure import *
from test_zams import *

test_chem()
print("chemistry passed testing")

test_adiabat()
print("adiabatic EOS passed testing")

test_pp()
print("reactions passed testing")

test_central_thermal()
print("central temperature passed testing")

test_Teff()
print("effective temperature passed testing")

test_surface_luminosity()
print("surface luminosity passed testing")
