"""
Code to plot the observed and modeled harmonic constituent
amplitude and phase.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun

pth = os.path.abspath('../../LiveOcean/plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import matplotlib.pyplot as plt
import pickle

import obsfun as ofn

home = os.environ.get('HOME')
dir00 = home + '/Documents/'
dir0 = dir00 + 'ptools_output/tide/'

noaa_sn_dict, dfo_sn_dict, sn_dict = ofn.get_sn_dicts()

# load observational data
year  = 2017

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
    
mod_dir = dir0 + 'mod_data/cas3_v0_lo6m/'
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
    
# plotting

plt.close('all')

# plot amplitude and phase, comparing two stations (ocean and inland)
# for both observations and model

# OFFSHORE
#name0 = 'Westport'
name0 = 'La Push'
#name0 = 'Neah Bay'
#name0 = 'Bamfield'
#name0 = 'Garibaldi'

# INLAND
#name1 = 'Campbell River'
#name1 = 'Point Atkinson'
#name1 = 'Vancouver'
name1 = 'Seattle'
#name1 = 'Seattle'
#name1 = 'Friday Harbor'

fig = plt.figure(figsize=(14, 8))
fig.suptitle('OCEAN:' + name0 + ' to INLAND:' + name1)

flo = .5
fhi = 2.5

ax1 = fig.add_subplot(221)
ax1.set_xlim(flo, fhi)
ax1.set_ylim(0, 3)
ax1.grid()
ax1.text(.05,.9, 'Amplification Factor (INLAND/OCEAN)', weight='bold', color='k',
    transform=ax1.transAxes)
ax1.text(.05,.8, 'OBSERVATION', weight='bold', color='r',
    transform=ax1.transAxes)
ax1.text(.05,.7, 'MODEL', weight='bold', color='b',
    transform=ax1.transAxes)

ax2 = fig.add_subplot(222)
ax2.set_xlim(flo, fhi)
ax2.set_ylim(0, 180)
ax2.grid()
ax2.text(.05,.9, 'Phase Shift (INLAND-OCEAN deg)', weight='bold', color='k',
    transform=ax2.transAxes)

hn_list = ['M2','S2','N2','O1','P1','K1']
for hn in hn_list:
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
    
ax3 = fig.add_subplot(223)
ax3.set_xlim(flo, fhi)
ax3.set_ylim(0, 1.3)
ax3.grid()
ax3.set_xlabel('Frequency (cycles/day)')
ax3.text(.05,.9, 'OCEAN Amplitude (m)', weight='bold', color='k',
    transform=ax3.transAxes)

ax4 = fig.add_subplot(224)
ax4.set_xlim(flo, fhi)
ax4.set_ylim(0, 1.3)
ax4.grid()
ax4.set_xlabel('Frequency (cycles/day)')
ax4.text(.05,.9, 'INLAND Amplitude (m)', weight='bold', color='k',
    transform=ax4.transAxes)

# the ocean phase looks good, so no need to plot
#ax4 = fig.add_subplot(224)
#ax4.set_xlim(flo, fhi)
#ax4.set_ylim(0, 360)
#ax4.grid()
#ax4.set_xlabel('Frequency (cycles/day)')
#ax4.text(.05,.9, 'OCEAN Phase (deg G)', weight='bold', color='k',
#    transform=ax4.transAxes)

hn_list = ['M2','S2','N2','O1','P1','K1']
for hn in hn_list:
    Ao0, Am0, Go0, Gm0, Fo0, Fm0 = get_AG(name0, hn, Hobs, Hmod)
    Ao1, Am1, Go1, Gm1, Fo1, Fm1 = get_AG(name1, hn, Hobs, Hmod)
    ax3.text(Fo0, Ao0, hn, color='r', weight='bold')
    ax3.text(Fm0, Am0, hn, color='b', weight='bold')
    ax4.text(Fo1, Ao1, hn, color='r', weight='bold')
    ax4.text(Fm1, Am1, hn, color='b', weight='bold')
    
if True:
    # Plot station locations
    fig2 = plt.figure(figsize=(8,12))
    ax = fig2.add_subplot(111)
    pfun.add_coast(ax, color='g')
    ax.set_xlim(-129, -121)
    ax.set_ylim(43, 51)
    pfun.dar(ax)
    for name in sn_dict.keys():
        xo = float(Mobs[name]['lon'])
        yo = float(Mobs[name]['lat'])
        ax.plot(xo, yo, '*r')
        ax.text(xo+.05, yo, name)
        xm = Mmod[name]['lon']
        ym = Mmod[name]['lat']
        ax.plot(xm, ym, '*b')
    
plt.show()