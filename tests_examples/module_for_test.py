"""
Test of access to modules from calling code.

NOTE: Testing this in Canopy can be tricky because when you first start
Canopy it runs ipython, which imports a bunch of things, including sys,
that are not present if you just start python from the command line. If you
then (in Canopy) do reset -f is will get rid of sys.  These concerns only
apply to the command line, not when invoking the code module_test, which
does not import extraneous stuff.

"""

import os
import sys

def add_to_path():   
    pth = os.path.abspath('../../LiveOcean/alpha')
    if pth not in sys.path:
        sys.path.append(pth)
    