"""
Code to experiment with different ROMS vertical coordinate parameters.
"""

# setup
import os; import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart(gridname='cas4', tag='v2')
# NOTE hmin = 4 in the cas4 grid

import zrfun
from importlib import reload
reload(zrfun)

import numpy as np
import matplotlib.pyplot as plt

# get S for our standard choices
s0 = Lfun.csv_to_dict(Ldir['grid'] + 'S_COORDINATE_INFO.csv')
S0 = zrfun.get_S(s0)
N = S0['N']

# original
# s0 = {'ITEMS': 'VALUES',
#  'THETA_S': '4',
#  'THETA_B': '0.8',
#  'TCLINE': '0',
#  'N': '30',
#  'VTRANSFORM': '1',
#  'VSTRETCHING': '1'}

# make alternate choices
s = s0.copy()
s['THETA_S'] = 4 # [0-10]
s['THETA_B'] = 2 # [0-4]
s['TCLINE'] = 10
s['VTRANSFORM'] = 2
s['VSTRETCHING'] = 4
S = zrfun.get_S(s)


NX = 100
x = np.arange(NX)
y = np.arange(N)
h = np.linspace(200, 5, NX)
z = zrfun.get_z(h, 0*h, S, only_w=True)

dz_deep = np.diff(z[:,0])
dz_shallow = np.diff(z[:,-1])

plt.close('all')
fig = plt.figure(figsize=(13,7))

ax = fig.add_subplot(231)
for rr in range(z.shape[0]):
    ax.plot(x,z[rr,:],'-b')
ax.set_title('Full Depth')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-h[0], 0)
ax.set_ylabel('Z (m)')

ax = fig.add_subplot(234)
for rr in range(z.shape[0]):
    ax.plot(x,z[rr,:],'-b')
ax.set_title('Near Surface')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-25, 0)
ax.set_ylabel('Z (m)')
ax.set_xlabel('X')
    
ax = fig.add_subplot(132)
ax.plot(dz_deep, y, '-*b')
ax.set_xlim(left=0)
ax.set_ylim(0,N-1)
ax.set_title('Depth = ' + str(h[0]) + ' (m)')
ax.set_xlabel('DZ (m) Deep')
ax.set_ylabel('Layer Number')
ii = 0
for item in ['theta_s', 'theta_b', 'tcline', 'N', 'Vtransform', 'Vstretching', 'hc']:
    ax.text(.05, 0.5 - .07*ii, item + ': ' + str(S[item]), transform=ax.transAxes, fontweight='bold')
    ii += 1
    
ax = fig.add_subplot(133)
ax.plot(dz_shallow, y, '-*b')
ax.set_xlim(left=0)
ax.set_ylim(0,N-1)
ax.set_title('Depth = ' + str(h[-1]) + ' (m)')
ax.set_xlabel('DZ (m) Shallow')
ax.set_yticklabels([])

plt.show()
