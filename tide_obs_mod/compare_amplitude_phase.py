"""
Compare observed and modeled amplitude and phase for all constituents
at all stations.  The goal is to determine optimum factors to use
when modifying the tidal forcing in the model.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import pickle
import numpy as np

import obsfun as ofn

dir0 = Ldir['parent'] + 'ptools_output/tide/'

# select model run
gtagex = 'cas4_v2_lo6biom'
year  = 2017
noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()
sn_list = sn_dict.keys()
    
def get_AG(hn, Hobs, Hmod):
    ho = Hobs
    hm = Hmod
    Ao = ho.A[ho.name==hn]
    Am = hm.A[hm.name==hn]
    Go = ho.g[ho.name==hn]
    Gm = hm.g[hm.name==hn]
    Fo = 24*ho.aux.frq[ho.name==hn] # cycles per day
    Fm = 24*hm.aux.frq[hm.name==hn]
    #
    return Ao, Am, Go, Gm, Fo, Fm

hn_list = ['M2','S2','N2','O1','P1','K1']

sn_coast = ['Charleston', 'South Beach', 'Garibaldi', 'Toke Point',
    'Westport', 'La Push', 'Neah Bay', 'Tofino', 'Bamfield']
sn_salish = ['Port Angeles', 'Friday Harbor', 'Cherry Point', 'Port Townsend',
    'Seattle', 'Tacoma', 'Point Atkinson', 'Vancouver', 'Patricia Bay',
    'Victoria Harbour', 'Campbell River', 'New Westminster']
    
df_dict = dict() # each DataFrame has one constituent
df = pd.DataFrame(index=sn_list, columns=['ar', 'dph'])
for hn in hn_list:
    df_dict[hn] = df.copy()

for name in sn_list:
    # load observational data
    obs_dir = dir0 + 'obs_data/'
    sn = sn_dict[name]
    hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    Hobs = pickle.load(open(hfn, 'rb'))
        
    # load model data
    mod_dir = dir0 + 'mod_data/' + gtagex + '/'
    sn = sn_dict[name]
    hfn = mod_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    Hmod = pickle.load(open(hfn, 'rb'))
    
    # get constituent info
    for hn in hn_list:
        Ao, Am, Go, Gm, Fo, Fm = get_AG(hn, Hobs, Hmod)
        df_dict[hn].loc[name, 'ar'] = Am/Ao
        # fix phase difference when they straddle 360
        if (Gm - Go) > 180:
            Gm = Gm - 360
        elif (Gm - Go) < -180:
            Gm = Gm + 360
        else:
            pass
        df_dict[hn].loc[name, 'dph'] = Gm - Go

print('\nCoast Stations: mean (std)')
for hn in hn_list:
    df = df_dict[hn]
    dff = df.loc[sn_coast,:]
    print(' %s: Amplitude Ratio = %5.2f (%5.5f), Phase Difference = %5.1f (%5.1f) [deg]' % (hn,
        dff['ar'].mean(), dff['ar'].std(), dff['dph'].mean(), dff['dph'].std()))

print('\nSalish Stations: mean (std)')
for hn in hn_list:
    df = df_dict[hn]
    dff = df.loc[sn_salish,:]
    print(' %s: Amplitude Ratio = %5.2f (%5.5f), Phase Difference = %5.1f (%5.1f) [deg]' % (hn,
        dff['ar'].mean(), dff['ar'].std(), dff['dph'].mean(), dff['dph'].std()))
