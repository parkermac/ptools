"""
Module to be run by module_test.py.

NOTE we are able to import modules at the start
that are then available to all functions defined in the module.

"""

import os
import sys

a = 10

def add_to_path():   
    pth = os.path.abspath('../../LiveOcean/alpha')
    if pth not in sys.path:
        sys.path.append(pth)

def change_dict(in_dict):
    in_dict['a'] = 'changed'
    
def change_var_noreturn(x):
    x = 2*x
    print(' - x in module = ' + str(x))
    
def change_var_return(x):
    x = 2*x
    print(' - x in module = ' + str(x))
    return x
    
def look_for_var():
    return a