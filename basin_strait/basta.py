"""
Code solve a system of basins and straits.
"""
import numpy as np

import bst_fun as bsf
from importlib import reload
reload(bsf)

pr = bsf.get_pr()

nd = 1000
NT_ic = int((nd+1)/(pr['dt']/86400))
NT = int((nd+1)/(pr['dt']/86400))

# specifiy which variable to use (always include 's')
v_list = ['s','dye']

# get the initial setup
Ba, St = bsf.make_BaSt('base', v_list)

# create output holders
Ba_out, Ba_out_cols, St_out, St_out_cols = bsf.make_BaSt_out(Ba, St, NT, v_list)

# run to steady state
fcase = 'base'
junk = bsf.bcalc(Ba, St, Ba_out, Ba_out_cols, St_out, St_out_cols, NT_ic, fcase, ic=True)

# run for real
Ba, St, Ba_out, St_out, TD = bsf.bcalc(Ba, St, Ba_out, Ba_out_cols,
        St_out, St_out_cols, NT, fcase)
# all will be time series except Ba and St, which are the final state

# plotting
import matplotlib.pyplot as plt
plt.close('all')
fig, axes = plt.subplots(3, 1, squeeze=False, figsize=(12,7))

ax = axes[0,0]
hatch_list = [ '/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*' ]
hh = 0
for bk in Ba.keys():
    if bk == 'a':
        pass
    else:
        S1 = Ba_out[:,Ba_out_cols[(bk,'S1')]]
        S2 = Ba_out[:,Ba_out_cols[(bk,'S2')]]
        ax.fill_between(TD, S1, S2, alpha=.5, label=bk, hatch=hatch_list[hh])
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
    ax.plot(TD, -St_out[:,St_out_cols[(sk,'q2')]]/1e3, '-', linewidth=3, label=sk)
    ax.grid(True)
    ax.set_xlim(TD[0], TD[-1])
aa = ax.get_ylim()
ax.set_ylim(aa[0]-1, aa[1]+1)
ax.legend(loc='upper right')
#ax.grid(True)
ax.set_xticklabels('')
ax.text(.05, .85, 'Strait Bottom Layer Q [$1000\ m^{3}s^{-1}$] (positive in)',
    fontsize=15, transform=ax.transAxes)
    
ax = axes[2,0]
for bk in Ba.keys():
    if bk == 'a':
        pass
    else:
        Qr = Ba_out[:,Ba_out_cols[(bk,'Qr')]]
        ax.plot(TD, Qr/1e3, '-', linewidth=3, label=bk)
        ax.set_xlim(TD[0], TD[-1])
ax.set_ylim(bottom=0)
ax.legend(loc='upper right')
#ax.grid(True)
ax.text(.05, .85, 'Basin River Flow [$1000\ m^{3}s^{-1}$]',
    fontsize=15, transform=ax.transAxes)
ax.set_xlabel('Time (days)')

plt.show()
