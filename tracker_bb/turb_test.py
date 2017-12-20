"""
Test to determine functionality of vertical turbulence. 
The result should not be a high concentration in the low diffusivity areas.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# number of particles
pnum = 100
# 1 hour timestep
delta_t = 3600
# number of days to run
dnum = 100

# create depths and AKs
zAKs = np.arange(-100,0,1)
AKs = -0.0005*np.cos(4*np.pi*-zAKs/100)+0.0005

# initial positions
pcs0 = -np.random.rand(pnum)
# change in pcs for a total of a 2m difference (water column is 100m)
dpcs = -0.01

pcs_df = pd.DataFrame(index=range(dnum*24), columns=range(pnum))
pcs_df.ix[0] = pcs0

# run hourly interpolation for 1 day
for j in np.arange(1,dnum*24):
    # set pcs
    pcs = np.array(pcs_df.ix[j-1], dtype='float64')
    
    # upper variables
    pcsu = pcs-dpcs
    pcsu[pcsu > 0] = 0
    AKsu = np.interp(pcsu, zAKs/100, AKs)
    
    # lower variables
    pcsb = pcs+dpcs
    pcsb[pcsb < -1] = -1
    AKsb = np.interp(pcsb,zAKs/100, AKs)
    
    # gradient between points 2m apart
    dAKs = (AKsu-AKsb)/2
    
    # update position with 1/2 of dAKs
    pcsm = pcs + (delta_t * dAKs/2)/100
    pcsm[pcsm > 0] = 0
    pcsm[pcsm < -1] = -1
    
    # updated AKs
    AKs1 = np.interp(pcsm, zAKs/100, AKs)
    
    # turbulence calculation from Banas, MacCready, and Hickey (2009)
    # w_turbulence = rand*sqrt(2k/delta(t)) + delta(k)/delta(z)
    # rand = random array with normal distribution
    rand = np.random.standard_normal(pnum)
    # only using half of gradient, first half already added
    turb = rand*np.sqrt(2*AKs1/delta_t) + dAKs/2
    
    # update position with rest of turbulence
    pcsf = pcsm + (delta_t * turb)/100
    pcsf[pcsf > 0] = 0
    pcsf[pcsf < -1] = -1
    
    # store final position
    pcs_df.ix[j] = pcsf


# Plotting

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

ax1.plot(AKs, zAKs, 'k')
ax1.set_ylim([-100,0])
ax1.set_xticklabels(['0', '2', '4', '6', '8', '10'])
ax1.set_title('AKs')
ax1.set_xlabel('Diffusivity ($10^{-4}  m^2s^{-1}$)')
ax1.set_ylabel('Depth (m)')


for j in pcs_df.columns:
    ax2.plot(np.arange(dnum*24)/24, pcs_df[j]*100, 'k', alpha=0.25, zorder=0)
ax2.scatter(np.zeros(pnum), pcs_df.iloc[0]*100, s=50, c='r', zorder=1)
ax2.scatter(np.ones(pnum)*dnum, pcs_df.iloc[-1]*100, s=50, c='r', zorder=1)
ax2.set_xlim([-dnum/25,dnum+dnum/25])
#ax2.set_xticklabels(np.arange(0,dnum*24,dnum*24/6))
ax2.set_title('Particles Over ' + str(dnum) + ' Days')
ax2.set_xlabel('Time (days)')


plt.show(fig)
