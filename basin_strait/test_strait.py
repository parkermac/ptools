"""
Code to test the strait function.
"""
import numpy as np
import matplotlib.pyplot as plt

import bst_fun as bsf
from importlib import reload
reload(bsf)

pr = bsf.get_pr()
v_list = ['s'] # set variables to use (always include 's')
case = 'base' # choose basin and strait setup
N = 100

# get the initial setup
Ba, St = bsf.make_BaSt(case, v_list)
sinfo = St['ba']

# set central point for parameter space exploration
Qr_0 = 1000
s1m_0 = 25
Ut_0 = 1

def reset_sinfo(sinfo, v_list, Qr_0, Ut_0, s1m_0):
    # get central state
    sinfo['Ut'] = Ut_0
    sinfo['s1m'] = s1m_0
    sinfo['s1p'] = pr['Socn']
    sinfo['s2m'] = pr['Socn']
    sinfo['s2p'] = pr['Socn']
    sinfo['q1'] = 1000
    sinfo['q2'] = -500
    sinfo['rev'] = False
    sinfo = bsf.strait(Qr_0, sinfo, v_list)
    return sinfo

def reset_s_vecs(N):
    s1m = np.nan * np.ones(N)
    s1p = np.nan * np.ones(N)
    s2m = np.nan * np.ones(N)
    s2p = np.nan * np.ones(N)
    q1 = np.nan * np.ones(N)
    q2 = np.nan * np.ones(N)
    return s1m, s1p, s2m, s2p, q1, q2
    
def store_s(sinfo, ii, s1m, s1p, s2m, s2p, q1, q2):
    if sinfo['rev'] == True:
        print('%d: reversed' % (ii))
    s1m[ii] = sinfo['s1m']
    s1p[ii] = sinfo['s1p']
    s2m[ii] = sinfo['s2m']
    s2p[ii] = sinfo['s2p']
    q1[ii] = sinfo['q1']
    q2[ii] = sinfo['q2']
    return s1m, s1p, s2m, s2p, q1, q2
    
# plotting, and experiments
plt.close('all')
fig, axes = plt.subplots(2, 3, squeeze=False, figsize=(12,7))
def plot_result(col, s1m, s1p, s2m, s2p, q1, q2, ivar, ivar_name):
    ax = axes[0,col]
    ax.plot(ivar, s1m, label='s1m')
    ax.plot(ivar, s1p, label='s1p')
    ax.plot(ivar, s2m, label='s2m')
    ax.plot(ivar, s2p, label='s2p')
    ax.set_xlim(ivar[0], ivar[-1])
    ax.set_ylim(0, 31)
    #ax.set_xlabel(ivar_name)
    if col==0:
        ax.set_ylabel('Strait End Salinities')
    ax.legend()
    
    # add lines defining the central values we used
    if ivar_name=='Forcing Salinity Gradient: s2p-s1m':
        ds = pr['Socn'] - s1m_0
        ax.plot([ds, ds], [0,31],'-k',linewidth=3,alpha=.3)
    elif ivar_name=='Qr (1000 m3/s)':
        ax.plot([Qr_0/1000, Qr_0/1000], [0,31],'-k',linewidth=3,alpha=.3)
    elif ivar_name=='Ut (m/s)':
        ax.plot([Ut_0, Ut_0], [0,31],'-k',linewidth=3,alpha=.3)
    
    ax = axes[1,col]
    ax.plot(ivar, q1/1e3, label='q1')
    ax.plot(ivar, -q2/1e3, label='-q2')
    ax.set_xlim(ivar[0], ivar[-1])
    ax.set_ylim(0, 25)
    ax.set_xlabel(ivar_name)
    if col==0:
        ax.set_ylabel('Layer transports (1000 m3/s)')
    ax.legend()
    

# EXPERIMENTS
# vary the salinity gradient
s1m, s1p, s2m, s2p, q1, q2 = reset_s_vecs(N)
s1m = np.linspace(0,pr['Socn'],N)
ii = 0
for this_s1m in s1m:
    sinfo = reset_sinfo(sinfo, v_list, Qr_0, Ut_0, this_s1m)
    sinfo = bsf.strait(Qr_0, sinfo, v_list)
    s1m, s1p, s2m, s2p, q1, q2 = store_s(sinfo, ii, s1m, s1p, s2m, s2p, q1, q2)
    ii += 1
DS = s2p - s1m
ivar = DS # set independent variable
ivar_name = 'Forcing Salinity Gradient: s2p-s1m'
plot_result(0, s1m, s1p, s2m, s2p, q1, q2, ivar, ivar_name)

# vary the river flow
s1m, s1p, s2m, s2p, q1, q2 = reset_s_vecs(N)
Qr = np.linspace(1,10e3,N)
ii = 0
for this_Qr in Qr:
    sinfo = reset_sinfo(sinfo, v_list, this_Qr, Ut_0, s1m_0)
    sinfo = bsf.strait(this_Qr, sinfo, v_list)
    s1m, s1p, s2m, s2p, q1, q2 = store_s(sinfo, ii, s1m, s1p, s2m, s2p, q1, q2)
    ii += 1
ivar = Qr/1000 # set independent variable
ivar_name = 'Qr (1000 m3/s)'
plot_result(1, s1m, s1p, s2m, s2p, q1, q2, ivar, ivar_name)

# vary the tidal forcing
s1m, s1p, s2m, s2p, q1, q2 = reset_s_vecs(N)
Ut = np.linspace(.1,10,N)
ii = 0
for this_Ut in Ut:
    sinfo = reset_sinfo(sinfo, v_list, Qr_0, this_Ut, s1m_0)
    sinfo['Ut'] = this_Ut
    sinfo = bsf.strait(Qr_0, sinfo, v_list)
    s1m, s1p, s2m, s2p, q1, q2 = store_s(sinfo, ii, s1m, s1p, s2m, s2p, q1, q2)
    ii += 1
ivar = Ut # set independent variable
ivar_name = 'Ut (m/s)'
plot_result(2, s1m, s1p, s2m, s2p, q1, q2, ivar, ivar_name)

plt.show()

