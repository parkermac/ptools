# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 09:10:17 2016

@author: PM5

This is the last step in the process, where you create and
populate LiveOcean_data/grids/[gridname]/.

"""

import shutil

from importlib import reload
import gfun; reload(gfun)
G = gfun.gstart()
import Lfun
Ldir = Lfun.Lstart(gridname=G['gridname'])

#%% select grid file

fn = gfun.select_file(G)
in_fn = G['gdir'] + fn

#%% copy to LiveOcean_data/grids

# make sure directories exist
Lfun.make_dir(Ldir['data'] + 'grids/', clean=False)
out_dir = Ldir['grid']
Lfun.make_dir(out_dir, clean=True)

# copy files
shutil.copyfile(in_fn, out_dir + 'grid.nc')
shutil.copyfile(G['gdir'] + 'river_info.csv', out_dir + 'river_info.csv')

#%% and put in the S coordinate info

# copy file
shutil.copyfile(Ldir['data'] + 'grids/S_COORDINATE_INFO_1.csv',
                out_dir + 'S_COORDINATE_INFO.csv')

# create the S.mat file
import subprocess
func = ( "Z_make_S(\'" + out_dir + "\')" )
cmd = Ldir['which_matlab']
run_cmd = [cmd, "-nodisplay", "-r", func, "&"]
# I dropped the "-nojvm" because it was throwing an error.  Maybe
# this is because I am using a newer matlab release.
cwd = Ldir['LO']+'shared/Z_functions/'
# I had to add the cwd= parameter in order to run in that directory.
proc = subprocess.Popen(run_cmd, cwd=cwd,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate()
# "out" is the screen output of the matlab code, and
# "err" stores arror messages from subprocess (I think).
#print(out.decode()) # will print "out" nicely
