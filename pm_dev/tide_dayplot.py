"""
Code to plot tides from NOAA csv files.

"""

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

dir00 = '/Users/PM5/Documents/'

# red in tide date
indir0 = dir00 + 'Classes/2017_Coastal_320/tide/'
fn = 'CO-OPS__9447130__hr.csv'
# Date Time, Water Level, Sigma, I, L
# 2016-01-01 00:00,1.302,0.000,0,0
# 2016-01-01 01:00,1.401,0.009,0,0
# and I think the time is UTC
df = pd.read_csv(indir0+fn, index_col='Date Time', parse_dates = True)
for k in df.keys():
    df = df.rename(columns={k: k.strip()})
df = df.drop(['Sigma', 'I', 'L'], axis=1)
df = df.rename(columns={'Water Level': 'SSH_obs'})

def make_dir(dirname, clean=False):
    # Make a directory if it does not exist.
    # Use clean=True to clobber the existing directory.
    import os
    if clean == True:
        import shutil
        shutil.rmtree(dirname, ignore_errors=True)
        os.mkdir(dirname)
    else:
        try:
            os.mkdir(dirname)
        except OSError:
            pass
            # assume OSError was raised because
            #directory already exists
            
# prepare a directory for results
outdir0 = dir00 + 'ptools_output/tide/'
make_dir(outdir0, clean=False)
outdir = outdir0 + 'dayplot/'
make_dir(outdir, clean=True)
    
# plotting loop

# initializing
plt.close('all')
d0 = datetime(2016,6,10,8,0,0) # start
# starting on hour 8 of UTC is midnight PST
d1 = d0 + timedelta(days=1)
d2 = d1 + timedelta(days=1)
# we get data from two days in order to do the blending
d_last = d0 + timedelta(days=(29.57 + 2)) # end
# 29.57 days is two S-N cycles
# but adding two days seems to get the ending phase back to
# the beginning phase - although the amplitudes don't match perfectly
# I assume because there are several other important constituents.
i_plot = 0

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1)
ax.grid()
ax.set_xlabel('Hour of Day (PST)')
ax.set_ylabel('Seattle Tide Height (m)')

while d2 <= d_last:    
    a = df.loc[d0:d1].values
    b = df.loc[d1:d2].values
    aa = np.array(a).flatten()
    bb = np.array(b).flatten()
    tt = np.arange(0,len(aa))
    ndiv = 6
    i_fr = 0
    while i_fr < ndiv:
        nouts = ('0000' + str(i_plot))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
        if (i_plot/10).is_integer():
            print(outname)
        fr = i_fr/ndiv
        cc = (1-fr)*aa + fr*bb # blended signal
        lh = ax.plot(tt, cc, '-b', linewidth=5)
        ax.set_title(d0.strftime('%Y.%m.%d'))
        ax.set_xlim(0, 24)
        ax.set_ylim(-1, 5)
        plt.draw()
        plt.savefig(outfile)
        # remove the line
        lh0 = lh.pop(0)
        lh0.remove()
        i_fr += 1
        i_plot += 1   
    d0 = d1
    d1 = d0 + timedelta(days=1)
    d2 = d1 + timedelta(days=1)

# and make a movie
import subprocess
cmd = ['ffmpeg',
    '-r', '24',
    '-i', outdir+'plot_%04d.png',
    '-vcodec', 'libx264',
    '-pix_fmt', 'yuv420p',
    '-crf', '25',
    outdir + 'movie.mp4']
proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if False:
    print('\n-main: screen output from subprocess-')
    print(proc.stdout.decode())
    print('\n-main: errors from subprocess-')
    # for some reason the ffmpeg output ends up in stderr
    print(proc.stderr.decode())

