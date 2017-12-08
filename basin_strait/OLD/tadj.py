"""
Adjustment time for the basin strait model.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import basin_strait_functions as bsf
from importlib import reload
reload(bsf)

dims, params = bsf.get_dims()
Socn = params['Socn']

ND = 3000 # number of days

# result vectors
vkm3_list = [20, 50, 100, 200, 500, 1000]

lb_list = [True, False]

for landward_basin in lb_list:
    
    # DataFrame for results
    tadj = pd.DataFrame(index=vkm3_list,
        columns=['S1', 'S2', 'Qin', 'v1', 'v2', 'qr0', 'qr1', 'qin0', 'qin1'])
        # the first three are adjustment times in days,
        # the rest are volumes, start/end river flow, start/end Qin
    
    for vkm3 in vkm3_list:
        dims['V1'] = vkm3*1e9 * 0.2 # volume of upper basin box
        dims['V2'] = vkm3*1e9 - dims['V1'] # volume of deeper basin box
        dims['V3'] = 0.2*vkm3*1e9 * 0.2 # volume of upper basin box
        dims['V4'] = 0.2*vkm3*1e9 - dims['V3'] # volume of deeper basin box
        # prepare result vectors
        # time
        dt = 86400 # time step (s)
        NT = int(ND*dt / 86400)
        T = np.linspace(0, NT*dt, NT)
        Td = T/86400
        # basin salinities
        S1 = np.nan * np.ones(NT); S2 = np.nan * np.ones(NT)
        S3 = np.nan * np.ones(NT); S4 = np.nan * np.ones(NT)
        # landward sill
        Qin = np.nan * np.ones(NT); Sin = np.nan * np.ones(NT)
        # seaward sill
        Qinp = np.nan * np.ones(NT); Sinp = np.nan * np.ones(NT)
        # forcing
        QR = 2000 * np.ones(NT)
        QR[0] = 1000
        UT = np.ones(NT)
        # draft intial conditions
        s1 = Socn-2; s2 = Socn-1; s3 = Socn-3; s4 = Socn-2
        # make an equilibrated initial condition
        for nt in range(NT-1):
            s1, s2, s3, s4, Qin[nt], Qinp[nt] = bsf.advance_s(
                    QR[0], UT[0], dims, params, s1, s2, s3, s4, dt,
                    landward_basin=landward_basin)
        # reset the initial condition to equilibrated state
        S1[0] = s1; S2[0] = s2; S3[0] = s3; S4[0] = s4
        # actual time integration
        for nt in range(NT-1):
            S1[nt+1], S2[nt+1], S3[nt+1], S4[nt+1], Qin[nt], Qinp[nt] = bsf.advance_s(
                    QR[nt+1], UT[nt+1], dims, params, S1[nt], S2[nt], S3[nt], S4[nt], dt,
                    landward_basin=landward_basin)
        # finishing up (we already have all the S1-4 final points)
        nt = NT-1
        junk, junk, junk, junk, Qin[nt], Qinp[nt] = bsf.advance_s(
                QR[nt], UT[nt], dims, params, S1[nt], S2[nt], S3[nt], S4[nt], dt,
                landward_basin=landward_basin)

        # function to calculate adjustment time
        def calc_tadj(name, ser, Td):
            v0 = ser[0]
            v1 = ser[-1]
            dv = v1-v0
            vv = (ser - v0)/dv
            value = (1 - 1/np.e)
            idx_q = (np.abs(vv-value)).argmin()
            idx_t0 = (np.abs(Td)).argmin()
            t0 = Td[idx_t0]
            t1 = Td[idx_q]
            tadj = t1-t0
            return tadj

        tadj.loc[vkm3,'v1'] = dims['V1']
        tadj.loc[vkm3,'v2'] = dims['V2']
        tadj.loc[vkm3,'qr0'] = QR[0]
        tadj.loc[vkm3,'qr1'] = QR[-1]
        tadj.loc[vkm3, 'qin0'] = Qin[0]
        tadj.loc[vkm3, 'qin1'] = Qin[-1]
    
        vdict = {'S1':S1, 'S2':S2, 'Qin':Qin}
        for vn in vdict.keys():
            this_tadj = calc_tadj(vn, vdict[vn], Td)
            tadj.loc[vkm3,vn] = this_tadj
            #print('%s: Adjustment Time = %0.2f days' % (vn, tadj))

    # calculate residence times
    tadj['tres1'] = (1/86400) * tadj['v1']/(tadj['qr0'] + tadj['qin0'])
    tadj['tres2'] = (1/86400) * tadj['v2']/tadj['qin0']

    if landward_basin == False:
        tadj_1b = tadj.copy()
    else:
        tadj_2b = tadj.copy()
    

# PLOTTING
fs_big = 16
fs_small = 12
lw_big = 3
lw_small = 2
bsf.set_rc(fs_big, fs_small, lw_big, lw_small)

#plt.close('all')
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

for landward_basin in lb_list:
    
    if landward_basin == False:
        tadj = tadj_1b.copy()
        alpha = 1
        color = 'c'
        tag = ''
        lw = 6
    else:
        tadj = tadj_2b.copy()
        alpha = 1
        color = 'r'
        tag = ' Two Basins'
        lw = 3
    mks = 20
    leg = True
    tadj.plot(y='S2', ax=ax, style='-', color=color, lw=lw,
        legend=leg, alpha=alpha, markersize=mks, label='Adjustment Time'+tag)
    tadj.plot(y='tres2', ax=ax, style=':', color=color, lw=lw,
        legend=leg, alpha=alpha, markersize=mks, label='Residence Time'+tag)
    
ax.set_xlabel('$Basin\ Volume\ (km^{3})$')
ax.set_ylabel('Days')
ax.grid()

ax.set_title('Deep Layer of Seaward Basin')

plt.show()

# RC CLEANUP
plt.rcdefaults()


    

