"""
Code to plot tides from NOAA csv files.

It makes a MOVIE where each frame is the tide height over a day,
so you can see how low tide marches forward avery day.

"""
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

home = os.environ.get('HOME')
dir00 = home + '/Documents/'

# prepare a directory for results
outdir0 = dir00 + 'ptools_output/tide/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + 'day_movie/'
Lfun.make_dir(outdir, clean=True)

# read in tide data (time is UTC)
indir = dir00 + 'ptools_data/tide/'
fn = 'CO-OPS__9447130__hr.csv'
# Date Time, Water Level, Sigma, I, L
# 2016-01-01 00:00,1.302,0.000,0,0
df = pd.read_csv(indir+fn, index_col='Date Time', parse_dates = True)
for k in df.keys():
    df = df.rename(columns={k: k.strip()})
df = df.drop(['Sigma', 'I', 'L'], axis=1)
df = df.rename(columns={'Water Level': 'SSH_obs'})

# PLOTTING LOOP

# initializing
plt.close('all')
d0 = datetime(2016,6,10,8,0,0) # start
# starting on hour 8 of UTC is midnight PST
d1 = d0 + timedelta(days=1)
d2 = d1 + timedelta(days=1)
# we get data from two days in order to do the blending

# set the number of days to plot
if False:
    # testing
    d_last = d0 + timedelta(days=5)
else:
    # 29.57 days is two S-N cycles
    d_last = d0 + timedelta(days=29)

i_plot = 0 # plot name number

# SET RC
fs1 = 16
fs2 = 20
lw1 = 3
lw2 = 5
plt.rc('xtick', labelsize=fs1)
plt.rc('ytick', labelsize=fs1)
plt.rc('xtick.major', size=10, pad=5, width=lw1)
plt.rc('ytick.major', size=10, pad=5, width=lw1)
plt.rc('axes', lw=lw1)
plt.rc('lines', lw=lw2)
plt.rc('font', size=fs2)
plt.rc('grid', color='g', ls='-', lw=lw1, alpha=.3)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1)
ax.grid()
ax.set_xlabel('Hour of Day (PST)')
ax.set_ylabel('Seattle Tide Height (m)')
ax.set_xlim(0, 24)
ax.set_ylim(-1, 5)

first = True # part of a scheme to make the looping smoother
while d1 <= d_last:
    a = df.loc[d0:d1].values
    if first:
        a0 = a.copy()
        first = False
    if d1 == d_last:
        b = a0
    else:
        b = df.loc[d1:d2].values
    aa = np.array(a).flatten()
    bb = np.array(b).flatten()
    tt = np.arange(0,len(aa))
    ndiv = 8 # number of divisions of the day
        # we do this to make a smoother movie
    i_fr = 0 # index into the divisions of the day
    while i_fr < ndiv:
        nouts = ('0000' + str(i_plot))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
        if (i_plot/10).is_integer():
            print(outname)
        fr = i_fr/ndiv
        cc = (1-fr)*aa + fr*bb # blended signal
        lh = ax.plot(tt, cc, '-b')
        ax.set_title(d0.strftime('%m/%d/%Y'))
        plt.tight_layout()
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

# RC CLEANUP
plt.rcdefaults()

# and make a movie
import subprocess
cmd = ['ffmpeg',
    '-r', '24', # framerate fps
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

