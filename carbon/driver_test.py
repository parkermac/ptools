"""
Code to test using CO2SYS.m

2020.07.08 Now including a test of the python version.

- Read in a roms history file.
- Save selected fields as matlab .mat files.
- Run a matlab subprocess to process these files.
- Read the results back into python and plot.

RESULTS (cas6 on my mac):

Field preparation takes 7 sec for a full layer, and 0.1 sec for Willapa region.

CO2SYS takes _FOREVER_ for one full layer, and 5 sec for just the Willapa region.

PyCO2SYS is much faster 10-36 sec for a full layer, and 1 sec for the Willapa region,
and the results were basically identical when I checked for the Willapa region

"""

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
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

from time import time
import cmocean.cm as cm

# (1) Get and package input fields
# file to work on
fn = (Ldir['roms'] + '/output/' + 'cas6_v3_lo8b/f2019.07.04/ocean_his_0001.nc')

do_both = False # only use the matlab version for a limited region
# specify a geographic region (aa=[] for full region)
if False:
    aa = [-124.4, -123.6, 46, 47.2] # Willapa Bay OA2 plot
else:
    aa = []
if len(aa) == 4:
    do_both = True

# specify layer number (S-coordinates, 0=bottom, -1=top)
NZ = 20
# get the fields
tt0 = time()
v_dict, plon, plat = cfun.get_layer(fn, NZ=NZ, aa=aa, print_info=True)
print('\n TIME to do calculation:')
print(' -- prepare fields %0.1f s' % (time()-tt0))

if do_both:
    # run CO2SYS
    tempdir = '../../ptools_output/carbon/temp/'
    PH, OM = cfun.get_carbon(v_dict, tempdir, print_info=True)

# run PyCO2SYS
tt0 = time()
from PyCO2SYS import CO2SYS
CO2dict = CO2SYS(v_dict['alkalinity'], v_dict['TIC'], 1, 2, v_dict['salt'], v_dict['temp'], v_dict['temp'],
    v_dict['pres'], v_dict['pres'], 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
# compare the above to this call from worker.m:
# A = CO2SYS(a.alkalinity(:), a.TIC(:), 1, 2, a.salt(:), a.temp(:), a.temp(:), a.pres(:), a.pres(:), 50, 2, 1, 10, 1);
#
# NOTE that the "in" and "out" versions of the returned variables will be identical because we pass the same
# pressure and temperature for input and output (in-situ in both cases)
print(' -- run PyCO2SYS %0.1f s' % (time()-tt0))
PH_alt = CO2dict['pHout']
PH_alt = PH_alt.reshape((v_dict['salt'].shape))
OM_alt = CO2dict['OmegaARout']
OM_alt = OM_alt.reshape((v_dict['salt'].shape))

if do_both:
    vdiff = .01 # limit for differences

    # PLOTTING
    plt.close('all')
    fig = plt.figure(figsize=(16,11))

    ax = fig.add_subplot(231)
    cs = plt.pcolormesh(plon, plat, PH[1:-1, 1:-1], vmin=7, vmax=8.5, cmap='Spectral')
    plt.colorbar(cs)
    ax.set_title('PH')
    pfun.dar(ax)

    ax = fig.add_subplot(232)
    cs = plt.pcolormesh(plon, plat, PH_alt[1:-1, 1:-1], vmin=7, vmax=8.5, cmap='Spectral')
    plt.colorbar(cs)
    ax.set_title('PH_alt')
    pfun.dar(ax)

    ax = fig.add_subplot(233)
    cs = plt.pcolormesh(plon, plat, PH_alt[1:-1, 1:-1]-PH_alt[1:-1, 1:-1], vmin=-vdiff, vmax=vdiff, cmap=cm.balance)
    plt.colorbar(cs)
    ax.set_title('Difference')
    pfun.dar(ax)

    ax = fig.add_subplot(234)
    cs = plt.pcolormesh(plon, plat, OM[1:-1, 1:-1], vmin=0, vmax=3, cmap='coolwarm_r')
    plt.colorbar(cs)
    ax.set_title('Omega')
    pfun.dar(ax)

    ax = fig.add_subplot(235)
    cs = plt.pcolormesh(plon, plat, OM_alt[1:-1, 1:-1], vmin=0, vmax=3, cmap='coolwarm_r')
    plt.colorbar(cs)
    ax.set_title('Omega_alt')
    pfun.dar(ax)

    ax = fig.add_subplot(236)
    cs = plt.pcolormesh(plon, plat, OM_alt[1:-1, 1:-1]-OM_alt[1:-1, 1:-1], vmin=-vdiff, vmax=vdiff, cmap=cm.balance)
    plt.colorbar(cs)
    ax.set_title('Difference')
    pfun.dar(ax)
    
else:
    
    # PLOTTING
    plt.close('all')
    fig = plt.figure(figsize=(16,11))

    ax = fig.add_subplot(121)
    cs = plt.pcolormesh(plon, plat, PH_alt[1:-1, 1:-1], vmin=7, vmax=8.5, cmap='Spectral')
    plt.colorbar(cs)
    ax.set_title('PH_alt')
    pfun.dar(ax)

    ax = fig.add_subplot(122)
    cs = plt.pcolormesh(plon, plat, OM_alt[1:-1, 1:-1], vmin=0, vmax=3, cmap='coolwarm_r')
    plt.colorbar(cs)
    ax.set_title('Omega_alt')
    pfun.dar(ax)

if 'mac' in Ldir['lo_env']:
    plt.show()
else:
    out_fn = tempdir + 'testplot.png'
    print('\nSAVING Figure to:')
    print(' -- ' + out_fn)
    plt.savefig(out_fn)
