"""
Code to look at NetCDf files interactively using imshow().

"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

# *************** USER INPUTS *******************************

# input file (ROMS history file)
fn = ('/Users/PM5/Documents/LiveOcean_roms/output/' +
    'cascadia1_base_lobio1/f2015.09.19/ocean_his_0002.nc')
    
# variable name
vn = 'PH'
    
# ***********************************************************
    
# set up the axes
plt.close('all')
fig = plt.figure(figsize=(16,10)) # (16,10) is good for my laptop
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2) # map
ax2 = plt.subplot2grid((1,3), (0,2), colspan=1) # buttons

# initialize the data plot
ds = nc.Dataset(fn)
# currently the default behavior is to assume we have a 3D field
# packed (z, y, x) and that we will cycle through z (sigma) levels.
vv = ds[vn][0,:,::-1,:] # reverse y-axis to work with imshow().
N = vv.shape[0]
k = 0 # initialize the slice number
V = vv[k,:,:]
cs = ax1.imshow(V, interpolation='nearest')
ax1.set_title(vn + ' Layer k = ' + str(k))
fig.colorbar(cs, ax=ax1)
plt.show()

# create control buttons
# list is organized from bottom to top
blist = ['playBack', 'stepBack', 'stepForward', 'playForward', 'done']
# nicer names for display
Blist = ['Play Backward', 'Step Backward', 'Step Forward', 'Play Forward', 'Done']
NB = len(blist) # number of buttons
ybc = np.arange(NB+1) - .5
offset = 1e5 # kludgey way to distinguish buttons from data field
xbc = np.arange(2) + offset
XBC, YBC = np.meshgrid(xbc, ybc)
ZB = np.arange(1,NB+1).reshape(NB,1)
cmap2 = plt.get_cmap(name='viridis')
ax2.pcolormesh(XBC,YBC,ZB, vmin=0, vmax=NB, cmap = cmap2)
ax2.set_axis_off()

def addButtonLabel(ax, plon, plat, nb, lab, tcol='k'):
    # draw and label buttons
    pad = .1
    ax.add_patch(
        plt.Rectangle((plon[0]+pad,plat[nb]+pad),
                      np.diff(plon)[-1]-2*pad, np.diff(plat)[-1]-2*pad,
                      fill=True, facecolor=inactive_color,
                      edgecolor='w'))
    ax.text(plon.mean(),nb, lab, fontsize=15,
             horizontalalignment='center', verticalalignment='center',
             color=tcol)
             
def setTitle(ax, vn, k):
    ax.set_title(vn + ': Layer k = ' + str(k))

# label the buttons (numbered bottom to top, 0 to NB-1)
bdict = dict(zip(range(NB), blist))
Bdict = dict(zip(range(NB), Blist))
active_color = 'k'
inactive_color = 'w'
for bnum in bdict.keys():
    addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=active_color)
        
# allow user to push the buttons to control data displayed
flag_get_ginput = True # Make False to exit the ginput loop

while flag_get_ginput:
    # get ginput
    a = plt.ginput(n=1, timeout=-1)
    # returns a list of tuples - of length 1
    b = np.array(a)
    b = np.round(b).astype(int)
    if b.shape != (1,2):
        b = np.array([[-1000, -1000]])
    ix = b[0, 0]
    iy = b[0, 1]

    # this code deals with button input
    if (ix >= offset):
        # were are in the buttons
        nb = iy # button number
        if (bdict[nb]=='playBack'):
            for k in range(N-1,-1,-1):
                V = vv[k,:,:]
                cs.set_data(V)
                setTitle(ax1, vn, k)
                plt.pause(.001)
        elif (bdict[nb]=='stepBack'):
            k -= 1
            if k <= 0:
                k = N-1
            V = vv[k,:,:]
            cs.set_data(V)
            setTitle(ax1, vn, k)
        elif (bdict[nb]=='stepForward'):
            k += 1
            if k >= N:
                k = 0
            V = vv[k,:,:]
            cs.set_data(V)
            setTitle(ax1, vn, k)
        elif (bdict[nb]=='playForward'):
            for k in range(N):
                V = vv[k,:,:]
                cs.set_data(V)
                setTitle(ax1, vn, k)
                plt.pause(.001)
        elif (bdict[nb]=='done'):
            flag_get_ginput = False
            ax1.set_title('DONE')
            for bnum in bdict.keys():
                addButtonLabel(ax2, xbc, ybc, bnum, Bdict[bnum], tcol=inactive_color)
        else:
            pass


