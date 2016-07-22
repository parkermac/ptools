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

#%% and put in the S coordinate info

shutil.copyfile(Ldir['res'] + 'S_COORDINATE_INFO.csv',
                Ldir['run'] + 'S_COORDINATE_INFO.csv')

import subprocess
func = ( "make_S(\'" + Ldir['gridname'] + "\')" )
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func, "&"]
# I had to add the cwd= parameter in order to run in that directory.
proc = subprocess.Popen(run_cmd, cwd=Ldir['res'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate() # "out" is the screen output of the matlab code
#print(out.decode()) # this ends up as part of the make_forcing_main.py screen output
