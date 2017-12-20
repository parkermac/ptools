"""
mort.py takes the output dictionary from tracking experiments and applies
a random, daily mortality to the particles.

The folder and filename to apply mortality to will be chosen in the run.

Output: Dictionary with the same shape as input, but with NaNs after 
particles have stopped.
"""

import pickle
import numpy as np

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

# percentage of particles to kill each day
mort = 25

indir = Ldir['LOo'] + 'tracks/'

# choose the run directory
print('\n%s\n' % '** Choose file to apply mortality **')
d_list_raw = os.listdir(indir)
d_list = []
for d in d_list_raw:
    d_list.append(d)
Ndt = len(d_list)
for ndt in range(Ndt):
    print(str(ndt) + ': ' + d_list[ndt])
my_ndt = int(input('-- Input number -- '))
dirname = d_list[my_ndt] + '/'

# create the list of run files
m_list_raw = os.listdir(indir + dirname)
m_list = []
for m in m_list_raw:
    if 'dailymort' not in m:
        m_list.append(m)
Npt = len(m_list)
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_ndt = int(input('-- Input number (99 for all) -- '))
if my_ndt == 99:
    pass
else:
    m_list = [m_list[my_ndt],]

# loop through each file
for inname in m_list:

    # compile list of day files
    p_list = os.listdir(indir + dirname + inname)
    p_list.sort()
    
    # loop through each day
    for p in p_list:
        
        # load data
        if p == p_list[0]:
            # day 0 contains P, Ldir, and grid components
            P, G, S, Ldir = pickle.load( open( indir + dirname + inname + '/' + p, 'rb' ) )
            NP = P['lon'].shape[1]
            # dictionary for location of dead particles
            mort_loc = dict()
            mort_loc['lon'] = np.zeros(NP)
            mort_loc['lat'] = np.zeros(NP)
            # array showing all particles that are still alive
            alive = np.ones(NP, dtype=bool)
        else:
            # other days only contain P and Ldir
            P, Ldir = pickle.load( open( indir + dirname + inname + '/' + p, 'rb' ) )
            
        # set P for particles killed on previous days
        # doesn't change 'ot', sets future locations to mort_loc, and other variables to nan
        dead = alive == False
        for var in P:
            if var in ['ot', 'age']:
                pass
            elif var in ['lon', 'lat']:
                P[var][:,dead] = mort_loc[var][dead]
            else:
                P[var][:,dead] = np.nan

        # only include already released larvae
        rel = P['age'] > 0
        
        # create random mask of all non-dead particles
        # note - because int rounds down, there will never be 0 particles left
        if p == p_list[0]:
            mort_mask = []
        else:
            mort_mask = np.random.choice(np.arange(NP)[alive*rel], int(sum(alive*rel)*(mort/100)))
        
        # update alive array
        alive[mort_mask] = False
        
        # store location at death
        mort_loc['lon'][mort_mask] = P['lon'][-1,mort_mask]
        mort_loc['lat'][mort_mask] = P['lat'][-1,mort_mask]
        

        # save the results
        outname = (inname + '_dailymort' + str(mort) + '/')
        outdir = indir + dirname + outname
        if p == p_list[0]:
            Lfun.make_dir(outdir, clean=True)
            pickle.dump((P, G, S, Ldir), open(outdir + p, 'wb'))
        else:
            pickle.dump((P, Ldir), open(outdir + p, 'wb'))
        print('Finished ' + p)

    print(50*'*')     
    print('Results saved to:\n' + outdir)
    print(50*'*')


