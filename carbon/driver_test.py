"""
Code to test using CO2SYS.m

- Read in a roms history file.
- Save selected fields as matlab .mat files.
- Run a matlab subprocess to process these files.
- Read the results back into python and plot.

RESULT: CO2SYS takes 14 sec to run for one full layer of the cas4 grid, on my mac.
It takes 8 sec for just the Willapa region. Almost all the time is taken by CO2SYS.m.
All other tasks take minimal time.
"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

if 'mac' in Ldir['lo_env']: # mac version
    pass
else: # regular (remote, linux) version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

import carbon_fun as cfun
from importlib import reload
reload(cfun)

# (1) Get and package input fields
# file to work on
fn = (Ldir['roms'] + '/output/' + 'cas4_v2_lo6biom/f2018.09.29/ocean_his_0001.nc')
# specify a geographic region (aa=[] for full region)
aa = [-124.4, -123.6, 46, 47.2] # Willapa Bay OA2 plot
# specify layer number (S-coordinates, 0=bottom, -1=top)
NZ = 0
# get the fields
v_dict, plon, plat = cfun.get_layer(fn, NZ=0, aa=aa, print_info=True)

# (2) run CO2SYS
tempdir = '../../ptools_output/carbon/temp/'
PH, OM = cfun.get_carbon(v_dict, tempdir, print_info=True)

# PLOTTING
plt.close('all')
fig = plt.figure(figsize=(16,8))

ax = fig.add_subplot(121)
cs = plt.pcolormesh(plon, plat, PH[1:-1, 1:-1], vmin=7, vmax=8.5, cmap='Spectral')
plt.colorbar(cs)
ax.set_title('PH')
pfun.dar(ax)

ax = fig.add_subplot(122)
cs = plt.pcolormesh(plon, plat, OM[1:-1, 1:-1], vmin=0, vmax=3, cmap='coolwarm_r')
plt.colorbar(cs)
ax.set_title('Omega')
pfun.dar(ax)

if 'mac' in Ldir['lo_env']:
    plt.show()
else:
    out_fn = tempdir + 'testplot.png'
    print('\nSAVING Figure to:')
    print(' -- ' + out_fn)
    plt.savefig(out_fn)
