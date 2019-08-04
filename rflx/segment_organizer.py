"""
This code gathers and organizes information about all the channel segments
in a realistic TEF/Efflux-Reflux (TEFx?) calculation.
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pickle

sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6',tag='v3')

sys.path.append(os.path.abspath(Ldir['LO'] + 'plotting'))
import pfun

sys.path.append(os.path.abspath(Ldir['parent'] + 'ptools/pgrid'))
import gfun
import gfun_plotting as gfp
Gr = gfun.gstart('cas6')

sys.path.append(os.path.abspath(Ldir['LO'] + 'x_tef'))
import tef_fun

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

indir0 = Ldir['LOo'] + 'tef/'
# choose the tef extraction to plot
if True:
    # hardwired choice for now
    item = 'cas6_v3_lo8b_2017.01.01_2017.12.31'
else:
    item = Lfun.choose_item(indir0)
indir = indir0 + item + '/thalweg/'
ThalMean = pickle.load(open(indir + 'ThalMean.p', 'rb'))
# ThalMean is the pickled results of x_tef process_thalweg_mean.py
# which is a dict consisting of the bulk 2-layer calculations
# ThalMean[ch_str] = (sect_list, q1, q2, qs1, qs2, s1, s2, dist)
# 1 = higer salinity, 2 = lower salinity


# segment definitions, a dict with entries consisting of
# Segment Name: ([list of surrounding sections], [list of rivers])
segs = {'J1':(['jdf1','jdf2'],['sanjuan', 'hoko'])}

