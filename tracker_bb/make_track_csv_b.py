# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:43:43 2016

@author: PM5
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import pandas as pd
import pickle  # python 3

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# choose the type of plot to make
# track directory
print('\n%s\n' % '** Choose mooring file to convert to csv **')
d_list = os.listdir(indir)
dNpt = len(d_list)
d_dict = dict(zip(range(dNpt), d_list))
for dnpt in range(dNpt):
    print(str(dnpt) + ': ' + d_list[dnpt])
my_d = d_dict[int(input('-- Input number (99 for all) -- '))]
# track file
m_list_raw = os.listdir(indir + my_d)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number (99 for all) -- '))

inname_list = []
if my_npt == 99:
    for ii in m_dict.keys():
        inname_list.append(m_dict[ii])   
else:
    inname_list.append(m_dict[my_npt])
    
for inname in inname_list:
      
    P, G, S, PLdir = pickle.load( open( indir + my_d + '/' + inname, 'rb' ) )
    
    # make sure the output directory exists
    outdir = Ldir['LOo'] + 'tracks/'
    outname = inname[:-2]
    odir = outdir + my_d + '/' + outname + '/'
    Lfun.make_dir(odir, clean=True)
    
    for vn in P.keys():
        a = pd.DataFrame(P[vn], index=P['ot'])
        a.index.name = 'seconds_since_1_1_1970'
        a.to_csv(odir + vn + '.csv')
