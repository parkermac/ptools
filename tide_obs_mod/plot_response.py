"""
Code to plot the observed and modeled harmonic constituent
amplitude and phase, for apair of stations (typically one
on the coast and one in the Salish Sea).

"""

# ========== USER ==========================

# (1) Select station pair

# OFFSHORE
#name0 = 'Tofino'
#name0 = 'Westport'
#name0 = 'La Push'
name0 = 'Neah Bay'
#name0 = 'Bamfield'
#name0 = 'Garibaldi'
#name0 = 'Victoria Harbour'

# INLAND
#name1 = 'Campbell River'
#name1 = 'Point Atkinson'
#name1 = 'Victoria Harbour'
name1 = 'Vancouver'
#name1 = 'Seattle'
#name1 = 'Tacoma'
#name1 = 'Friday Harbor'

# (2) Select model run

#gtagex = 'cas4_v2_lo6biom'
#gtagex = 'cas5_v3_lo8'
gtagex = 'cas6_v2_lo8'

# (3) Select year
year  = 2017

# =========================================


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

import matplotlib.pyplot as plt
import pickle

import obsfun as ofn

dir00 = Ldir['parent']
dir0 = dir00 + 'ptools_output/tide/'

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()

obs_dir = dir0 + 'obs_data/'
Tobs = dict()
Mobs = dict()
Hobs = dict()
for name in sn_dict.keys():
    # observed tide
    sn = sn_dict[name]
    fn = obs_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = obs_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = obs_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    #Tobs[name] = pd.read_pickle(fn)
    Mobs[name] = Lfun.csv_to_dict(mfn)
    Hobs[name] = pickle.load(open(hfn, 'rb'))
    
mod_dir = dir0 + 'mod_data/' + gtagex + '/'
Tmod = dict()
Mmod = dict()
Hmod = dict()
for name in sn_dict.keys():
    # observed tide
    sn = sn_dict[name]
    fn = mod_dir + 'tide_' + str(sn) + '_' + str(year) + '.p'
    mfn = mod_dir + 'm_' + str(sn) + '_' + str(year) + '.csv'
    hfn = mod_dir + 'h_' + str(sn) + '_' + str(year) + '.p'
    #Tmod[name] = pd.read_pickle(fn)
    Mmod[name] = Lfun.csv_to_dict(mfn)
    Hmod[name] = pickle.load(open(hfn, 'rb'))
    
def get_AG(name, hn, Hobs, Hmod):
    ho = Hobs[name]
    hm = Hmod[name]
    Ao = ho.A[ho.name==hn]
    Am = hm.A[hm.name==hn]
    Go = ho.g[ho.name==hn]
    Gm = hm.g[hm.name==hn]
    Fo = 24*ho.aux.frq[ho.name==hn] # cycles per day
    Fm = 24*hm.aux.frq[hm.name==hn]
    #
    return Ao, Am, Go, Gm, Fo, Fm
    
# PLOTTING
#plt.close('all')
fig, axes = plt.subplots(nrows=3, ncols = 3, squeeze=False, figsize=(14,12))
fig.suptitle(gtagex + ' OCEAN:' + name0 + ' to INLAND:' + name1)

# frequecy limits (cpd)
flo = .5
fhi = 2.5

# fontsize and weight
fs = 6
fw = 'normal'

ax1 = axes[0,0]
ax1.set_xlim(flo, fhi)
ax1.set_ylim(0, 3)
ax1.grid()
ax1.text(.05,.9, 'Amplification Factor (INLAND/OCEAN)',
    color='k',
    transform=ax1.transAxes)
ax1.text(.05,.8, 'OBSERVATION',
    color='r',
    transform=ax1.transAxes)
ax1.text(.05,.7, 'MODEL',
    color='b',
    transform=ax1.transAxes)

ax2 = axes[0,1]
ax2.set_xlim(flo, fhi)
ax2.set_ylim(0, 180)
ax2.grid()
ax2.text(.05,.9, 'Phase Shift (INLAND-OCEAN deg)',
    color='k',
    transform=ax2.transAxes)

for hn in ofn.hn_list:
    Ao0, Am0, Go0, Gm0, Fo0, Fm0 = get_AG(name0, hn, Hobs, Hmod)
    Ao1, Am1, Go1, Gm1, Fo1, Fm1 = get_AG(name1, hn, Hobs, Hmod)
    Aro = Ao1/Ao0
    Arm = Am1/Am0
    dGo = Go1 - Go0
    if dGo < 0:
        dGo += 360
    dGm = Gm1 - Gm0
    if dGm < 0:
        dGm += 360
    ax1.plot(Fo0, Aro, '*r', Fm0, Arm, '*b')
    ax2.plot(Fo0, dGo, '*r', Fm0, dGm, '*b')

# AMPLITUDES
ax3 = axes[1,0]
ax3.set_xlim(flo, fhi)
ax3.set_ylim(0, 1.4)
ax3.grid()
#ax3.set_xlabel('Frequency (cycles/day)')
ax3.text(.05,.9, 'OCEAN Amplitude (m)',
    color='k',
    transform=ax3.transAxes)
#
ax4 = axes[2,0]
ax4.set_xlim(flo, fhi)
ax4.set_ylim(0, 1.4)
ax4.grid()
ax4.set_xlabel('Frequency (cycles/day)')
ax4.text(.05,.9, 'INLAND Amplitude (m)',
    color='k',
    transform=ax4.transAxes)
#
for hn in ofn.hn_list:
    Ao0, Am0, Go0, Gm0, Fo0, Fm0 = get_AG(name0, hn, Hobs, Hmod)
    Ao1, Am1, Go1, Gm1, Fo1, Fm1 = get_AG(name1, hn, Hobs, Hmod)
    ax3.text(Fo0, Ao0, hn, color='r', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    ax3.text(Fm0, Am0, hn, color='b', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    ax4.text(Fo1, Ao1, hn, color='r', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    ax4.text(Fm1, Am1, hn, color='b', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')

# PHASES
ax3 = axes[1,1]
ax3.set_xlim(flo, fhi)
ax3.set_ylim(0,360)
ax3.grid()
#ax3.set_xlabel('Frequency (cycles/day)')
ax3.text(.05,.9, 'OCEAN Phase (deg)',
    color='k',
    transform=ax3.transAxes)
#
ax4 = axes[2,1]
ax4.set_xlim(flo, fhi)
ax4.set_ylim(0,360)
ax4.grid()
ax4.set_xlabel('Frequency (cycles/day)')
ax4.text(.05,.9, 'INLAND Phase (deg)',
    color='k',
    transform=ax4.transAxes)
#
for hn in ofn.hn_list:
    Ao0, Am0, Go0, Gm0, Fo0, Fm0 = get_AG(name0, hn, Hobs, Hmod)
    Ao1, Am1, Go1, Gm1, Fo1, Fm1 = get_AG(name1, hn, Hobs, Hmod)
    ax3.text(Fo0, Go0, hn, color='r', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    ax3.text(Fm0, Gm0, hn, color='b', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    ax4.text(Fo1, Go1, hn, color='r', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    ax4.text(Fm1, Gm1, hn, color='b', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')

# AMPLITUDE RATIOS
ax5 = axes[1,2]
ax5.set_xlim(flo, fhi)
ax5.set_ylim(0.5, 1.7)
ax5.grid()
#ax5.set_xlabel('Frequency (cycles/day)')
ax5.text(.05,.9, 'OCEAN Amplitude Ratio (Mod/Obs)',
    color='k',
    transform=ax5.transAxes)
#
ax6 = axes[2,2]
ax6.set_xlim(flo, fhi)
ax6.set_ylim(0.5, 1.7)
ax6.grid()
ax6.set_xlabel('Frequency (cycles/day)')
ax6.text(.05,.9, 'INLAND Amplitude Ratio (Mod/Obs)',
    color='k',
    transform=ax6.transAxes)
#
for hn in ofn.hn_list:
    Ao0, Am0, Go0, Gm0, Fo0, Fm0 = get_AG(name0, hn, Hobs, Hmod)
    Ao1, Am1, Go1, Gm1, Fo1, Fm1 = get_AG(name1, hn, Hobs, Hmod)
    ax5.text(Fo0, Am0/Ao0, hn+' '+str(int(100*Am0/Ao0)/100),
            color='k', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    ax6.text(Fo1, Am1/Ao1, hn+' '+str(int(100*Am1/Ao1)/100),
            color='k', weight=fw, fontsize=fs, horizontalalignment='center', verticalalignment='center')
    
ax_map = axes[0,2]
# Plot station locations
pfun.add_coast(ax_map, color='g')
ax_map.set_xlim(-129, -121)
ax_map.set_ylim(46, 51)
pfun.dar(ax_map)
for name in sn_dict.keys():
    xo = float(Mobs[name]['lon'])
    yo = float(Mobs[name]['lat'])
    if name in [name0, name1]:
        ax_map.plot(xo, yo, '*r')
        ax_map.text(xo+.05, yo, name)
    else:
        ax_map.plot(xo, yo, 'ok', alpha=.3)
        
    
plt.show()