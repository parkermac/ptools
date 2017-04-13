# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:18:09 2016

@author: PM5

Organizational functions for pgrid.
"""

# USER EDIT

gridname = 'cas1'

dir0 = '/Users/PM5/Documents/'
pgdir = dir0 + 'ptools_output/pgrid/'

if 'aestus' in gridname:
    ri_dir = dir0 + 'ptools_output/river/analytical/'
else:
    ri_dir = dir0 + 'ptools_output/river/pnw_all_2016_07/'

# END USER EDIT

import os
import sys
alp = os.path.abspath(dir0 +'LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

plp = os.path.abspath(dir0 +'LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)

def gstart():
    gdir = pgdir + gridname + '/'
    Gr ={'gridname': gridname, 'dir0': dir0, 'pgdir': pgdir, 'gdir': gdir,
         'ri_dir': ri_dir}
    return Gr

def select_file(Gr, using_old_grid=False):
    # interactive selection
    if using_old_grid==True:
        fn_list = []
        dir0 = '/Users/PM5/Documents/LiveOcean_data/grids/'
        gn_list = ['cascadia1', 'cascadia2']
        for gn in gn_list:
            fn_list.append(dir0 + gn + '/grid.nc')
    elif using_old_grid==False:
        print('\n** %s in <<%s>> **\n' % ('Choose file to edit', Gr['gridname']))
        fn_list_raw = os.listdir(Gr['gdir'])
        fn_list = []
        for item in fn_list_raw:
            if item[-3:] == '.nc':
                fn_list.append(item)
    Nfn = len(fn_list)
    fn_dict = dict(zip(range(Nfn), fn_list))
    for nfn in range(Nfn):
        print(str(nfn) + ': ' + fn_list[nfn])
    my_nfn = int(input('-- Input number -- '))
    fn = fn_dict[my_nfn]
    return fn

def increment_filename(fn, tag='_m'):
    # create the new file name
    gni = fn.find(tag)
    new_num = ('00' + str(int(fn[gni+2: gni+4]) + 1))[-2:]
    fn_new = fn.replace(fn[gni:gni+4], tag + new_num)

    return fn_new

def default_choices(Gr):
    # Default choices (can override in each case)    
    dch = dict()    
    # Set analytical to true when we define the bathymetry analytically.
    dch['analytical'] = False    
    # Adjustment to zero of the bathymetry to account for the fact
    # that mean sea level is somewhat higher than NAVD88.
    dch['use_z_offset'] = True
    dch['z_offset'] = -1.06    
    # Set this to True for grids that are much coarser than the bathy files.
    # It means we average all data inside a rho grid cell to get depth there,
    # instead of using interpolation.
    # We also make a mask_rho that is the average of all values in a cell
    # using 1=water, 0=land.
    dch['do_cell_average'] = False    
    # With fjord_cliff_edges=True the smoothing deviates from its
    # usual volume-conserving
    # nature when it is next to a masked region, and instead adjusts the slope
    # by preferentially deepening at the coast.  This does a much better job of
    # preserving thalweg depth in channels like Hood Canal
    dch['fjord_cliff_edges'] = True    
    # Decide if the grid will allow wetting and drying
    dch['wet_dry'] = False    
    # Set the minimum depth, and decide if it should be enforced.
    dch['use_min_depth'] = False
    dch['min_depth'] = 4 # meters, positive down    
    # specify topography files to use
    dch['t_dir'] = Gr['dir0'] + 'tools_data/geo_data/topo/'
    # list of topo files: coarsest to finest
    dch['t_list'] = ['srtm15/topo15.grd',
              'cascadia/cascadia_gridded.mat',
             'psdem/PS_183m.mat',
             'ttp_patch/TTP_Regional_27m_patch.mat']
    return dch


    



