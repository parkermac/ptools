#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 09:55:11 2018

@author: pm7

Test of argparse.

RESULT: For my typical useage I think it is best to only define
the default and the type.  Then argparse will throw an error if:
    - a flag is given with no value
    - the value does not match the type
and if a flag is not given it will be filled with the default.

Note: you need to use boolean_string() to parse booleans correctly.

"""

import argparse

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

def prog(in_str):
    parser = argparse.ArgumentParser()

    # RESULT: For my typical useage I think it is best to just define
    # the default and the type.  Then argparse will throw an error if:
    #     - a flag is given with no value
    #     - the value does not match the type
    # and if a flag is not given it will be filled with the default.
    parser.add_argument('-s', '--string', default='hi', type=str)
    parser.add_argument('-i', '--integer', default=10, type=int)    

    # this way returns False with no flag, and True with the flag
    parser.add_argument('-ft', '--flag_true', action='store_true')

    # this way returns True with no flag, and False with the flag
    parser.add_argument('-ff', '--flag_false', action='store_false')    
    
    # to get a boolean this way is more complicated because is requires you
    # to tell it either True or False, and it requires the boolean_string function,
    # but it is clear and robust, returning a true boolean
    parser.add_argument('-b', '--bb', default=True, type=boolean_string)    
    
    args = parser.parse_args(in_str.split())
    d = args.__dict__     
    return d

in_list = ['-s hello', '-i 5', '-b False', '-ft', '-ff', '-b True', '-b False']
#in_list = ['-s a,b,c']

for in_str in in_list:
    print('\nInput string: ' + in_str)
    d = prog(in_str)
    print(d)