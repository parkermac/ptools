# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:18:09 2016

@author: PM5

Utility function for pgrid.
"""

gridname = 'test'

dir0 = '/Users/PM5/Documents/'
pgdir = dir0 + 'ptools_output/pgrid/'
gdir = pgdir + gridname + '/'

import os
import sys
# "alp" should also be on the path for any calling code
alp = os.path.abspath(dir0 +'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)

import matfun

def gstart():

    return gridname, dir0, pgdir, gdir

def increment_filename(fn, tag='_m'):
    # create the new file name
    gni = fn.find(tag)
    new_num = ('00' + str(int(fn[gni+2: gni+4]) + 1))[-2:]
    fn_new = fn.replace(fn[gni:gni+4],'_m' + new_num)

    return fn_new

def get_coast():
    c_dir = dir0 + 'tools_data/geo_data/coast/'
    c_file = 'pnw_coast_combined.mat'
    c_fn = c_dir + c_file
    cmat = matfun.loadmat(c_fn)

    return cmat



