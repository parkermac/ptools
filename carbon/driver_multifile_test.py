"""
Code to test using CO2SYS.m for a list of files.

RESULT for the Willapa region:

TIME to gather and process input fields:
 -- 73 input files
 -- 33.9 seconds

 TIME to do calculation:
 -- run CO2SYS 20.4 s

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import numpy as np
import time

import carbon_fun as cfun
from importlib import reload
reload(cfun)

# (1) Get and package input fields
# input dir to work on
indir = Ldir['roms'] + '/output/' + 'cas4_v2_lo6biom/f2018.09.29/'
hlist = os.listdir(indir)
hlist = [item for item in hlist if 'ocean_his' in item]
hlist.sort()
#hlist = hlist[:3]
fn_list = []
for hh in hlist:
    fn_list.append(indir + hh)
    
# specify a geographic region (aa=[] for full region)
aa = [-124.4, -123.6, 46, 47.2] # Willapa Bay OA2 plot
# specify layer number (S-coordinates, 0=bottom, -1=top)
NZ = 0

ii = 0
tt0 = time.time()
for fn in fn_list:
    # get the fields
    v_dict, plon, plat = cfun.get_layer(fn, NZ=0, aa=aa)
    if ii == 0:
        vv_dict = v_dict.copy()
        for vn in v_dict.keys():
            v = v_dict[vn]
            NR, NC = v.shape
            vv_dict[vn] = v.reshape((1, NR, NC))
    else:
        for vn in v_dict.keys():
            v = v_dict[vn]
            NR, NC = v.shape
            vv_dict[vn] = np.concatenate((vv_dict[vn], v.reshape((1, NR, NC))), axis=0)
    ii += 1
print('\nTIME to gather and process input fields:')
print(' -- %d input files' % (len(fn_list)))
print(' -- %0.1f seconds' % (time.time() - tt0))

# (2) run CO2SYS
tempdir = '../../ptools_output/carbon/temp/'
PH, OM = cfun.get_carbon(vv_dict, tempdir, print_info=True)

