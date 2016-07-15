# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 09:10:17 2016

@author: PM5

This is the last step in the process, where you send a selected grid file
and the river_info.csv file to LiveOcean/preamble/.

"""

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()

import os
import sys
fpth = os.path.abspath('../')
if fpth not in sys.path:
    sys.path.append(fpth)
import Lfun
Ldir = Lfun.Lstart(gridname=G['gridname'])

import shutil


#%% select grid file
fn = gfun.select_file(G)
in_fn = G['gdir'] + fn

#%% copy to LiveOcean

Lfun.make_dir(Ldir['run'], clean=False)

shutil.copyfile(in_fn, Ldir['run'] + 'grid.nc')
shutil.copyfile(G['gdir'] + 'river_info.csv', Ldir['run'] + 'river_info.csv')
