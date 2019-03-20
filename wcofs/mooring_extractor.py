"""
Extract multiple mooring-like records from WCOFS.
The input is structured to conform with layer_extractor.py.

"""

# setup
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
import numpy as np
from datetime import datetime, timedelta
start_time = datetime.now()
import zfun
import zrfun
import netCDF4 as nc

from importlib import reload

pth = os.path.abspath('../../LiveOcean/x_moor')
if pth not in sys.path:
    sys.path.append(pth)
import moor_fun as mfun
reload(mfun)

import wcofs_fun as wfun
reload(wfun)

# Load the module that defines mooring lists.
import moor_lists as ml
reload(ml)

job_name = 'comt3_2014_offshore'

# command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-g', '--gridname', nargs='?', type=str, default='WCOFS')
parser.add_argument('-t', '--tag', nargs='?', type=str, default='avg')
parser.add_argument('-x', '--ex_name', nargs='?', type=str, default='Exp29')
# see alpha/Lfun.get_fn_list() for acceptable list_type values

args = parser.parse_args()

# get list of history files to plot
Ldir = Lfun.Lstart(args.gridname, args.tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + args.ex_name

in_dir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'
fn_list_raw = os.listdir(in_dir)
fn_list = [(in_dir + ff) for ff in fn_list_raw if 'ocean_avg' in ff]
fn_list.sort()

# get info
fn = fn_list[0]
ds = nc.Dataset(fn)

G, S, T = wfun.get_basic_info(fn)
# for vn in ds.variables:
#     print(vn)

if False:
    zeta = ds['zeta'][:].squeeze()
    z_rho, z_w = zrfun.get_z(G['h'], zeta, S)

sta_dict, v2_list, v3_list_rho, v3_list_w = ml.get_sta_dict('comt3_2014_offshore')

sta = 'CA015'
xy = sta_dict[sta]
slon = xy[0]
slat = xy[1]

x, y = wfun.wcofs_lonlat_2_xy(slon, slat)

# plotting
import matplotlib.pyplot as plt
plt.close('all')

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.pcolormesh(G['x_rho'], G['y_rho'], G['h'])

plt.show()


