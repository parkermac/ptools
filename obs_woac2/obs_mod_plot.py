"""
This plots results of obs_mod_compare.py.

"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3', ex_name='lo8b')

# get the list of cruises and casts
sta_df_dir = Path(__file__).absolute().parent.parent.parent / 'ptools_output' / 'woac2'
OMa = pd.read_pickle(sta_df_dir / 'OMa.p')

vv_list = ['s', 'th', 'do', 'din', 'ta', 'dic', 'ph', 'arag']

plt.close('all')
fig = plt.figure(figsize=(14, 8))
ii = 1
for vn in vv_list:
    vnm = vn + '_m'
    ax = fig.add_subplot(2,4,ii)
    vmin = np.min((OMa[vn].min(),OMa[vnm].min()))
    vmax = np.max((OMa[vn].max(),OMa[vnm].max()))
    OMa.plot(x=vn, y=vnm, legend=False, xlim=(vmin,vmax), ylim=(vmin,vmax),
        ax=ax, style='.',markersize=2)
    ax.plot([vmin,vmax],[vmin,vmax],'-k')
    ax.axis('square')
    ii += 1
plt.show()

