"""
Code to process ORCA data files and save DataFrames.

Parker MacCready
"""
# setup
import os
import sys
path_to_alpha = '/Users/PM5/Documents/LiveOcean/alpha'
alp = os.path.abspath(path_to_alpha)
if alp not in sys.path:
    sys.path.append(alp)
from importlib import reload
import matfun
reload(matfun)

import numpy as np
import pandas as pd
from datetime import datetime, timedelta

# choose station(s) to plot
sta_name_list = [
    'HC_NB',
    'HC_DB',
    'HC_DK',
    'HC_HP',
    'HC_TW',
    'PW',
    'CI',
    ]

# choose variables to process
data_to_plot = [
    'Salinity',
    'Temperature',
    'Sigma',
    'Fluor',
    'DO',
    #'NO3'  # unreliable, don't use  
    ]

# station names available : and their long names
sta_name_dict = {'CI':'South Sound - Carr Inlet', 'HC_DB':'Hood Canal - Dabob Bay',
    'HC_DK':'Hood Canal - Duckabush', 'HC_HP':'Hood Canal - Hoodsport',
    'HC_NB':'Hood Canal - North Buoy', 'HC_TW':'Hood Canal - Twanoh',
    'PW':'Main Basin - Point Wells'}

# where are the data csv files
dir0 = '/Users/PM5/Documents/'
indir = dir0 + 'tools_data/obs_data/mooring/ORCA_new_2015_02/'
outdir0 = dir0 + 'ptools_output/env_data/'
outdir = outdir0 + 'ORCA/'
try:
    os.mkdir(outdir0)
except OSError:
    pass # assume OSError was raised because directory already exists
try:
    os.mkdir(outdir)
except OSError:
    pass # assume OSError was raised because directory already exists

# data names, units, and plot ranges
data_long_names=['Bsal',   'Btemp',       'Bsigma',     'Bfluor','Boxy_mgL',    'Bnitrate']
data_names =  ['Salinity', 'Temperature', 'Sigma',      'Fluor', 'DO',          'NO3']
data_units =  ['psu',      '$^{\circ}C$', '$kg m^{3}$', 'units', '$mg l^{-1}$', '$\mu M$']

# suppress warning messages (e.g. from YY>Ylo where YY is an array with nan's)
import warnings
warnings.simplefilter("ignore")

# dictionaries to look up data attributes using names
data_name_dict = dict(zip(data_long_names, data_names))
data_unit_dict = dict(zip(data_names, data_units))

#%%

for sta_name in sta_name_list:
       
    # load the data
    filename = sta_name + '_CTD_data_bin_web'
    fn = indir + filename + '.mat'
    a = matfun.loadmat(fn)
    A = dict()
    Z = a['Bpress'] # packed top to bottom
    ZZ = np.nanmean(Z, axis = 1)
    ZZind = np.arange(len(ZZ))
    T = a['Btime']
    tvec = np.nanmean(T, axis = 0)
    NZ, NT = Z.shape
    for long_name in data_name_dict:
        var_name = data_name_dict[long_name]
        A[var_name] = a[long_name]
        A[var_name][A[var_name]==-555.] = np.nan
    
    # make a time axis
    t0 = datetime(2000,1,1)
    TT = []
    for t in tvec:
        try:
            TT.append(t0 + timedelta(t))
        except ValueError:
            TT.append(np.nan)
    
    # intitialize data frames to save time series data
    dft = pd.DataFrame(index=TT, columns=data_to_plot)
    dfb = pd.DataFrame(index=TT, columns=data_to_plot)
   
    for dtp in data_to_plot:
        iz0 = ZZind[ZZ>=-5.].tolist()
        x0 = A[dtp][iz0,:]
        dft[dtp] = np.nanmean(x0, axis=0)
        
        iz1 = ZZind[ZZ<-5.].tolist()
        x1 = A[dtp][iz1,:]
        dfb[dtp] = np.nanmean(x1, axis=0)
    
    dft = dft.dropna()
    dft = dft.drop_duplicates() 
    dft = dft.sort()
    dft = dft[pd.notnull(dft.index)]
    
    dfb = dfb.dropna()
    dfb = dfb.drop_duplicates() 
    dfb = dfb.sort()
    dfb = dfb[pd.notnull(dfb.index)]
    
    # resample daily
    dftm = dft.resample('D', how='mean')
    dfbm = dfb.resample('D', how='mean')
    
    print('Saving ' + sta_name)
    sys.stdout.flush()
    dftm.to_pickle(outdir + 'df_top_' + sta_name + '.p')
    dfbm.to_pickle(outdir + 'df_bot_' + sta_name + '.p')




    


            

    

