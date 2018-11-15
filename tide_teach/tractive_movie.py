"""
Code to test coordinate rotation in 3D.

Notes on coordinates:

Longitude (-pi to pi) and latitude are just like on Earth,
i.e. longitude increases counter-clockwise.

The view angles azimuth and elevation are defined exactly
like lon and lat.

x,y,z are a right-handed coordinate system with:
x aligned with lon = 0
y aligned with lon = pi/2 (90 deg)
z aligned with the north pole

So when azim = 0 you are looking directly
from the direction that the x-axis is pointing
"""
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from importlib import reload
import tractive_functions as tf
reload(tf)

# SETUP

do_movie = True

# timekeeping vectors
mhour_vec = range(100) # this controls the number of frames
nmh = len(mhour_vec)
mlon_deg_vec = np.linspace(180, -180, nmh)
mlon_rad_vec = np.deg2rad(mlon_deg_vec)
#
r = 10 # radius of the Earth Sphere
mdec_deg = 23.5 # set the lunar declination
mdec_rad = np.deg2rad(mdec_deg)
#

if do_movie:
    home = os.environ.get('HOME')
    dir00 = home + '/Documents/'
    # prepare a directory for results
    outdir0 = dir00 + 'ptools_output/tide/'
    Lfun.make_dir(outdir0, clean=False)
    outdir = outdir0 + 'tractive_movie_' + str(mdec_deg) +'/'
    Lfun.make_dir(outdir, clean=True)

# PLOTTING

plt.close('all')

if do_movie == False:
    mhour_vec = [mhour_vec[int(nmh/4)]]
    
i_plot = 0    
for mh in mhour_vec:
    mlon_rad = mlon_rad_vec[mh]

    if do_movie:
        # name the output file
        nouts = ('0000' + str(i_plot))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
    
    fig = plt.figure(figsize=(9, 9)) # (9,9)
    ax = fig.add_subplot(111, projection='3d')
    
    tf.draw_sphere(ax, r)
    tf.draw_moon(ax, r, mlon_rad, mdec_rad)
    tf.draw_vector_ring(ax, mdec_deg, mlon_rad, r, direction='toward_moon')
    tf.draw_vector_ring(ax, mdec_deg, mlon_rad, r, direction='away_from_moon')
    tf.draw_tractive_pretzels(ax, mdec_rad, mlon_rad_vec, mlon_rad, r)
    #tf.add_axes(ax, r)
    
    # add_axes(ax)
    ax.set_axis_off()
    scl = 1#.92
    ax.set_xlim(-scl*r, scl*r)
    ax.set_ylim(-scl*r, scl*r)
    ax.set_zlim(-scl*r, scl*r)
    ax.set_aspect('equal')
    if do_movie:
        ax.azim = 20
        ax.elev = 14
    else:
        ax.azim = 0
        ax.elev = 0
    
    ax.set_title(('Lunar Declination = %s Degrees' % str(mdec_deg)),
        fontweight='bold', fontsize=16)
        
    fig.tight_layout()

    if do_movie:
        plt.draw()
        plt.savefig(outfile)
        plt.close()
    else:
        plt.show()
    
    i_plot += 1

if do_movie:
    tf.make_movie(outdir)
