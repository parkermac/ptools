"""
Code solve a system of basins and straits.
"""
import numpy as np

import basta_fun as bsf
from importlib import reload
reload(bsf)

Ba, St = bsf.get_BaSt('base')

# run to steady state
fcase = 'base'
junk = bsf.bcalc(St, Ba, 3001, fcase, ic=True)
# run for real
S11, S22, Qrr, Q11, Q22, Utt, TD, St, Ba = bsf.bcalc(St, Ba, 1001, fcase)

# plotting
import matplotlib.pyplot as plt
plt.close('all')
fig, axes = plt.subplots(2, 1, squeeze=False, figsize=(12,7))

ax = axes[0,0]
hatch_list = [ '/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*' ]
hh = 0
for bk in Ba.keys():
    if bk == 'a':
        pass
    else:
        S111 = S11[bk]
        S222 = S22[bk]
        ax.fill_between(TD, S111, S222, alpha=.5, label=bk, hatch=hatch_list[hh])
        ax.set_xlim(TD[0], TD[-1])
        hh += 1
aa = ax.get_ylim()
ax.set_ylim(aa[1]+1, aa[0]-1)
ax.legend(loc='upper right')
#ax.grid(True)
ax.set_xticklabels('')
ax.text(.05, .85, 'Basin Salinities',
    fontsize=15, transform=ax.transAxes)

ax = axes[1,0]
for sk in St.keys():
    ax.plot(TD, -Q22[sk]/1e3, '-', linewidth=3, label=sk)
    ax.grid(True)
    ax.set_xlim(TD[0], TD[-1])
aa = ax.get_ylim()
ax.set_ylim(aa[0]-1, aa[1]+1)
ax.legend(loc='upper right')
#ax.grid(True)
ax.text(.05, .85, 'Strait Bottom Layer Q [$1000\ m^{3}s^{-1}$] (positive in)',
    fontsize=15, transform=ax.transAxes)

plt.show()
