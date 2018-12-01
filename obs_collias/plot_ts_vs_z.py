"""
Plot Collias and Ecology CTD casts in Hood Canal to see if
s(z) really changed as much as my bulk analysis indicated.

"""

# imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# data locations
dir0 = '../../ptools_data/collias/'
dir00 = '../../ptools_data/ecology/'
# load processed station info and data
sta_df_ecology = pd.read_pickle(dir00 + 'sta_df.p')

dir1 = '../../ptools_output/collias/'

# create dicts that associate a region number (Collias system) with a name
# and with a list of modern Ecology station letters
collias_region_names = {1:'Strait of Juan de Fuca',
        2:'Admiralty Inlet',
        3:'Main Basin',
        4:'South Sound',
        5:'Hood Canal',
        6:'Whidbey Basin',
        7:'North Sound',
        8:'San Juan Island Passages'}
ecology_region_lists = {1:['SJF'],
        2:['ADM', 'PTH'],
        3:['PSB', 'ELB', 'EAP', 'CMB'],
        4:['OAK', 'BUD', 'CSE', 'DNA', 'CRR', 'GOR', 'NSQ'],
        5:['HCB'],
        6:['SKG', 'SAR', 'PSS'],
        7:['GRG', 'BLL'],
        8:['RSR']}
        
# ============================================================
region_list = range(1,9)#[3,5]
testing = False
yd0 = 120 # start of May
yd1 = 273 # Start of October
# ============================================================

# variables to process
col_list = ['Station', 'Date', 'Temp. (deg C)', 'Salinity', 'Z (m)']
plt.close('all')

for region in region_list:
    
    # First get raw Collias Data
    
    if testing:
        year_list = [1952]
    else:
        year_list = range(1932, 1976) # data go from 1932 to 1975
    
    print('\nCollias: working on region ' + str(region))
    Collias = pd.DataFrame()
    for year in year_list:
        print(' - year = ' + str(year))
        b = pd.read_pickle(dir0 + 'Bottles_' + str(year) + '.p')

        if len(b) > 0:
            b = b[b['Num0']==str(region)] # select region

            # now we want to go through all individual casts
            sta_list = b['Station'].values
            for sta in set(sta_list):
                bb = b[b['Station']==sta]
                date_list = bb['Date'].values
                for date in set(date_list):
                    c = bb[bb['Date']==date]
                    c = c.set_index('Z (m)')
                    c = c.reindex(columns=col_list)
                    # now we have a little DataFrame that is a single cast
                    c.loc[:,'Date'] = date
                    Collias = pd.concat((Collias,c))

    if len(Collias) > 0:
        Collias.loc[:,'Z (m)'] = Collias.index.values
        Collias = Collias.set_index('Date')

    # Next get raw Ecology Data
    if testing:
        year_list = [2006]
    else:
        year_list = range(1999, 2018)

    print('\nEcology: working on region ' + str(region))
    Ecology = pd.DataFrame() # DataFrame to hold all cast data
    
    for year in year_list:
        print(' - year = ' + str(year))
        Casts = pd.read_pickle(dir00 + 'Casts_' + str(year) + '.p')
        # loop over all stations
        sta_list = [sta for sta in sta_df_ecology.index if sta[:3] in ecology_region_lists[region]]
        for station in sta_list:
            casts = Casts[Casts['Station'] == station]
            casts = casts.set_index('Date')
            cast_vn_dict = {'Salinity': 'Salinity', 'Temperature': 'Temp. (deg C)', 'Z': 'Z (m)'}
            # rename the columns
            casts = casts.rename(columns=cast_vn_dict)
            # limit the columns
            casts = casts.reindex(columns=['Station', 'Z (m)', 'Salinity', 'Temp. (deg C)'])
            Ecology = pd.concat((Ecology,casts), sort=False)
            
    # Binning
    dz = 5
    z = np.arange(-150,dz,dz)
    zbins = z[:-1] + np.diff(z)/2

    EC_dry = pd.DataFrame(index=zbins, columns=['S Col', 'S Ecy', 'T Col', 'T Ecy'])
    for zb in zbins:
        a = Collias['Z (m)']>(zb-dz/2)
        b = Collias['Z (m)']<=(zb+dz/2)
        aa = Collias.index.dayofyear > yd0
        bb = Collias.index.dayofyear < yd1
        c = a & b & aa & bb
        EC_dry.loc[zb, 'S Col'] = Collias.loc[c, 'Salinity'].mean()
        EC_dry.loc[zb, 'T Col'] = Collias.loc[c, 'Temp. (deg C)'].mean()
    for zb in zbins:
        a = Ecology['Z (m)']>(zb-dz/2)
        b = Ecology['Z (m)']<=(zb+dz/2)
        aa = Ecology.index.dayofyear > yd0
        bb = Ecology.index.dayofyear < yd1
        c = a & b & aa & bb
        EC_dry.loc[zb, 'S Ecy'] = Ecology.loc[c, 'Salinity'].mean()
        EC_dry.loc[zb, 'T Ecy'] = Ecology.loc[c, 'Temp. (deg C)'].mean()
    EC_dry.index.name = 'Z (m)'
    EC_dry['Z (m)'] = EC_dry.index.values

    EC_wet = pd.DataFrame(index=zbins, columns=['S Col', 'S Ecy', 'T Col', 'T Ecy'])
    for zb in zbins:
        a = Collias['Z (m)']>(zb-dz/2)
        b = Collias['Z (m)']<=(zb+dz/2)
        aa = Collias.index.dayofyear < yd0
        bb = Collias.index.dayofyear > yd1
        c = a & b & (aa | bb)
        EC_wet.loc[zb, 'S Col'] = Collias.loc[c, 'Salinity'].mean()
        EC_wet.loc[zb, 'T Col'] = Collias.loc[c, 'Temp. (deg C)'].mean()
    for zb in zbins:
        a = Ecology['Z (m)']>(zb-dz/2)
        b = Ecology['Z (m)']<=(zb+dz/2)
        aa = Ecology.index.dayofyear < yd0
        bb = Ecology.index.dayofyear > yd1
        c = a & b & (aa | bb)
        EC_wet.loc[zb, 'S Ecy'] = Ecology.loc[c, 'Salinity'].mean()
        EC_wet.loc[zb, 'T Ecy'] = Ecology.loc[c, 'Temp. (deg C)'].mean()
    EC_wet.index.name = 'Z (m)'
    EC_wet['Z (m)'] = EC_wet.index.values

    # PLOTTING
    vn_list = ['Temp. (deg C)', 'Salinity', 'DO (mg L-1)', 'NO3 (uM)']
    lim_list = [[6,18],         [20,34],     [0,14],        [0,34]]
    lim_dict = dict(zip(vn_list, lim_list))
    
    ylim=(-150,0)

    fig = plt.figure(figsize=(12,8))
    
    # SALINITY
    xt = .05; yt = .8

    ax = fig.add_subplot(321)
    Collias.plot(x='Salinity',y='Z (m)',style='ob',markersize=2,alpha=.1, ax=ax, legend=False)
    Ecology.plot(x='Salinity',y='Z (m)',style='.r',markersize=.3,alpha=.1, ax=ax, legend=False)
    ax.set_xlim(lim_dict['Salinity'])
    ax.set_xlabel('')
    ax.set_xticklabels([])
    ax.set_ylabel('Z (m)')
    ax.grid(True)
    ax.set_ylim(ylim)
    ax.text(xt,yt,'(a) Raw Salinity', fontweight='bold', transform=ax.transAxes)

    ax = fig.add_subplot(323)
    EC_dry.plot(x='S Col',y='Z (m)',style='-ob',markersize=2,alpha=1, ax=ax, label='Collias: Dry Months')
    EC_dry.plot(x='S Ecy',y='Z (m)',style='-or',markersize=2,alpha=1, ax=ax, label='Ecology: Dry Months')
    ax.set_xlim(lim_dict['Salinity'])
    ax.set_xlabel('')
    ax.set_xticklabels([])
    ax.set_ylabel('Z (m)')
    ax.grid(True)
    ax.legend(loc='lower left')
    ax.set_ylim(ylim)
    ax.text(xt,yt,'(b) Salinity in 5 m z-bins', fontweight='bold', transform=ax.transAxes)

    ax = fig.add_subplot(325)
    EC_wet.plot(x='S Col',y='Z (m)',style='-ob',markersize=2,alpha=1, ax=ax, label='Collias: Wet Months')
    EC_wet.plot(x='S Ecy',y='Z (m)',style='-or',markersize=2,alpha=1, ax=ax, label='Ecology: Wet Months')
    ax.set_xlim(lim_dict['Salinity'])
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Z (m)')
    ax.grid(True)
    ax.legend(loc='lower left')
    ax.set_ylim(ylim)
    ax.text(xt,yt,'(c) Salinity in 5 m z-bins', fontweight='bold', transform=ax.transAxes)
    
    # TEMPERATURE
    xt = .95; yt=.8
    
    ax = fig.add_subplot(322)
    Collias.plot(x='Temp. (deg C)',y='Z (m)',style='ob',markersize=2,alpha=.1, ax=ax, legend=False)
    Ecology.plot(x='Temp. (deg C)',y='Z (m)',style='.r',markersize=.3,alpha=.1, ax=ax, legend=False)
    ax.set_xlim(lim_dict['Temp. (deg C)'])
    ax.set_xlabel('')
    ax.set_xticklabels([])
    ax.grid(True)
    ax.set_ylim(ylim)
    ax.text(xt,yt,'(d) Raw Temperature', fontweight='bold', transform=ax.transAxes, horizontalalignment='right')
    
    ax = fig.add_subplot(324)
    EC_dry.plot(x='T Col',y='Z (m)',style='-ob',markersize=2,alpha=1, ax=ax, label='Collias: Dry Months')
    EC_dry.plot(x='T Ecy',y='Z (m)',style='-or',markersize=2,alpha=1, ax=ax, label='Ecology: Dry Months')
    ax.set_xlim(lim_dict['Temp. (deg C)'])
    ax.set_xlabel('')
    ax.set_xticklabels([])
    ax.grid(True)
    ax.legend(loc='lower right')
    ax.set_ylim(ylim)
    ax.text(xt,yt,'(e) Temperature in 5 m z-bins', fontweight='bold', transform=ax.transAxes, horizontalalignment='right')

    ax = fig.add_subplot(326)
    EC_wet.plot(x='T Col',y='Z (m)',style='-ob',markersize=2,alpha=1, ax=ax, label='Collias: Wet Months')
    EC_wet.plot(x='T Ecy',y='Z (m)',style='-or',markersize=2,alpha=1, ax=ax, label='Ecology: Wet Months')
    ax.set_xlim(lim_dict['Temp. (deg C)'])
    ax.set_xlabel('Temp. (deg C)')
    ax.grid(True)
    ax.legend(loc='lower right')
    ax.set_ylim(ylim)
    ax.text(xt,yt,'(f) Temperature in 5 m z-bins', fontweight='bold', transform=ax.transAxes, horizontalalignment='right')
    
    fig.suptitle('Collias (1932-1975) & Ecology (1999-2017) Data: ' + collias_region_names[region])
    plt.savefig(dir1 + 'TSZ_' + collias_region_names[region].replace(' ','_') + '.png')


plt.show()

