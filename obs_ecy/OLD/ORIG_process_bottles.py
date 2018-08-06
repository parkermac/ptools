"""
Process data from WA Dept. of Ecology Bottle casts, and
save the results for future use.

"""

import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

# read in stored station info (links station names to locations)
dir0 = '../../ptools_data/ecology/'
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# read in bottle data
bottle_fn = dir0 + 'raw/Parker_2006-present_Nutrients.xlsx'
sheet_name = '2006-2017'
allbot = pd.read_excel(bottle_fn, sheet_name=sheet_name)

# select a single station
sbc_dict = dict()
for sn in sta_df.index:

    sba = allbot[allbot['Station']==sn]
    # set the index
    sba = sba.set_index('Date')

    # keep only selected columns
    sb = sba.reindex(columns=['NomDepth', 'Sampling Depth', 'PO4(uM)D', 'SiOH4(uM)D',
           'NO3(uM)D', 'NO2(uM)D', 'NH4(uM)D'])

    # deal with missing data
    sb[sb==999] = np.nan
    sb = sb.dropna()

    #sbn = sb.reindex(columns=['NO3(uM)D', 'NO2(uM)D', 'NH4(uM)D'])

    sb['DIN(uM)'] = sb.reindex(columns=['NO3(uM)D', 'NO2(uM)D', 'NH4(uM)D']).sum(axis=1)

    sb_dict = dict()

    # pull out only selected depth
    d_list = [0, 10, 30]
    for dd in d_list:
        sb_dict[dd] = sb[sb['NomDepth']==dd]

    # drop repeated values
    for dd in d_list:
        sb_dict[dd] = sb_dict[dd][~sb_dict[dd].index.duplicated()]

    # develop a monthly climatology
    sbc = pd.DataFrame(index=np.arange(1,13), columns=d_list)
    sbc.index.name = 'Month'

    for m in sbc.index:
        for dd in d_list:
            sbc.loc[m, dd] = sb_dict[dd].loc[sb_dict[dd].index.month==m, 'DIN(uM)'].mean()
    
    sbc_dict[sn] = sbc

# plot data
plt.close('all')
fig = plt.figure(figsize=(13,8))
ii = 0
NR = 5
NC = 8
for sn in sta_df.index:
    
    # print some infor for model initialization
    jmax = sbc_dict[sn].loc[1,30]
    des = sta_df.loc[sn,'Descrip']
    print('%s: January max DIN(uM) = %0.1f, (%s)' % (sn, jmax, des))
    ax = fig.add_subplot(NR,NC,ii+1)
    ir = int(np.floor(ii/NC))
    ic = int(ii - NC*ir)
    if sn == 'SJF000':
        do_leg = True
    else:
        do_leg = False
    sbc_dict[sn].plot(ax=ax, legend=do_leg, xticks=[3,6,9])
    if ic > 0:
        ax.set_ylabel('')
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('DIN(uM)')
    if ir < NR-1:
        ax.set_xlabel('')
        ax.set_xticklabels([])
    ax.set_xlim(1,12)
    ax.set_ylim(0,40)
    ax.grid(True)
    ax.text(.05,.85,sn, transform=ax.transAxes)
    d0 = des.split('-')[0].strip()
    ax.text(.05,.2,d0,fontsize=8, transform=ax.transAxes)
    try:
        d1 = des.split('-')[1].strip()
        ax.text(.05,.07,d1,fontsize=8, transform=ax.transAxes)
    except IndexError:
        pass
    ii += 1

# #vn_list = ['PO4(uM)D', 'SiOH4(uM)D','DIN(uM)', 'Sampling Depth']
# vn_list = ['DIN(uM)']
# NR = len(vn_list)
#
# rr = 1
# for vn in vn_list:
#     ax = fig.add_subplot(NR,1,rr)
#     for dd in d_list:
#         sb_dict[dd].plot.line(y=vn, style='-*', ax=ax, label=(str(dd) + ' m'))
#     ax.set_title(sn + ' ' + vn)
#     ax.grid(True)
#     ax.set_xlim(pd.datetime(2006,1,1), pd.datetime(2018,1,1))
#     if rr < NR:
#         ax.set_xlabel('')
#         ax.set_xticklabels([])
#     rr += 1
#
plt.show()
