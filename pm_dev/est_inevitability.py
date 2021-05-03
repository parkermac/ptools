"""
Code to explore the inevitability of estuarine structure given
the requirement that Ri = 1/4.
"""

import numpy as np
import matplotlib.pyplot as plt

import cubic_solver as cs
from importlib import reload
reload(cs)

Socn = 30
g = 9.8
beta = 7.7e-4
H = 20
Ric = 1/4

c = np.sqrt(g*beta*H*Socn)
Ubar = .01

bbar = Ubar/c

# for Ric = 1/4
if Ric == 1/4:
    print('using Ric=1/4 case')
    Da = cs.cubic_solver(1, -9*bbar**2, 36*bbar**2, -36*bbar**2)
    Db = 2*np.sqrt(Da)
else:
# general case for any Ric
    Da = cs.cubic_solver(1, -(36*Ric)*bbar**2, (144*Ric)*bbar**2, -(144*Ric)*bbar**2)
    Db = (1/np.sqrt(Ric))*np.sqrt(Da)

# translate to dimensional quantities
DS = Da*Socn
Sbar = Socn - DS/2
Du = Db * c

# make z-axis
N = 1001
z = np.linspace(0,H,N)
dz = np.diff(z)
zm = z[:-1] + dz/2 # z on midpoints (for integrals)
Nm = len(zm)

# generate u and z on zm
um = Du*zm/H - Du/2 + Ubar
sm = Socn - DS*zm/H

# check values
Ubar_check = um.mean()
Sbar_check = sm.mean()
F_check = 100 * (um*sm).mean() / (um.mean()*sm.mean())
Ri_check = g*beta*H*DS/(Du**2)
print('Ubar = %0.4f, Ubar_check = %0.4f' % (Ubar, Ubar_check))
print('Sbar = %0.4f, Sbar_check = %0.4f' % (Sbar, Sbar_check))
print('F = %0.2f , F_check = %0.2f%% of Ubar*Sbar' % (0, F_check))
print('Ri = %0.4f, Ri_check = %0.4f' % (Ric, Ri_check))

# plotting
plt.close('all')
fs = 14
plt.rc('font', size=fs)
fig = plt.figure(figsize=(14,7))

ax = fig.add_subplot(121)
ax.plot(um,zm, '-k')
ax.set_ylim(0,H)
ax.grid(True)
ax.set_xlabel('U (m/s)')
ax.set_ylabel('Z (m)')

ax = fig.add_subplot(122)
ax.plot(sm,zm, '-k')
ax.set_ylim(0,H)
ax.grid(True)
ax.set_xlabel('Salinity (g/kg)')

plt.show()
plt.rcdefaults()

# part 2, try adding an unmixed layer on the bottom

zzm = np.concatenate((zm-H, zm))
u0 = Ubar - Du/2
s0 = Socn
uum = np.concatenate((u0*np.ones(Nm), um))
ssm = np.concatenate((s0*np.ones(Nm), sm))

# adjust uum so that it satisfies  volume conservation
uum = uum - uum.mean() + Ubar

# adjust ssm so that it satisfies salt conservation
sadj = -(uum*ssm).mean() / uum.mean()
#ssm += sadj

# check values
Ubar_check = uum.mean()
F_check = 100 * (uum*ssm).mean() / (uum.mean()*ssm.mean())
print('')
print('Ubar = %0.4f, Ubar_check = %0.4f' % (Ubar, Ubar_check))
print('F = %0.2f , F_check = %0.2f%% of Ubar*Sbar' % (0, F_check))

# plotting
#plt.close('all')
fs = 14
plt.rc('font', size=fs)
fig = plt.figure(figsize=(14,7))

ax = fig.add_subplot(121)
ax.plot(uum,zzm, '-k')
ax.set_ylim(-H,H)
ax.grid(True)
ax.axvline(0)
ax.set_xlabel('U (m/s)')
ax.set_ylabel('Z (m)')

ax = fig.add_subplot(122)
ax.plot(ssm,zzm, '-k')
ax.set_ylim(-H,H)
ax.grid(True)
ax.set_xlabel('Salinity (g/kg)')

plt.show()
plt.rcdefaults()


