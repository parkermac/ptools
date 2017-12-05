# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 08:06:36 2016

@author: PM5

Code to automate getting multiple years of NDBC buoy data
"""

import bs4
import urllib.request as U

from importlib import reload
import ndbc_fun as ndf
reload(ndf)

# START USER EDITS

# Specify which NDBC buoy number(s) to get data for

get_all = True

if get_all:
    # all of them
    name_dict = ndf.get_name_dict()
    sn_list = [k for k in name_dict.keys()]
    yr_list = range(1970, 2020)
else:
    # or some of them (select by hand)
    #sn_list = ['46005', '46041', '46002', '46015', '46059', '46014']
    #sn_list = ['42040', '42035']
    sn_list = ['sisw1']
    yr_list = range(2015, 2018)

# END USER EDITS

# make sure ouput directory structure exists
data_dir00 = '../../ptools_data/'
data_dir0 = data_dir00 + 'ndbc/'
ndf.make_dir(data_dir00)
ndf.make_dir(data_dir0)

# get the data and save as text files for each station/year
for sn in sn_list:
    print('\n*** Getting ' + sn + ' ***\n')
    outdir = data_dir0 + sn + '/'
    ndf.make_dir(outdir, clean=True)
    

    for yr in yr_list:
        try:
            idn = str(sn) + 'h' + str(yr)
            url_str = ('http://www.ndbc.noaa.gov/view_text_file.php?filename=' +
                       idn + '.txt.gz&dir=data/historical/stdmet/')
            html = U.urlopen(url_str, timeout=10)
            soup = bs4.BeautifulSoup(html, 'html.parser')
            sn_text = soup.findAll(text=True)
            sns = str(sn_text)[2:-2]
            sns = sns.replace('\\n','\n')
            f = open(outdir + idn + '.txt','w')
            f.write(sns)
            f.close()
            print(' Retrieved ' + idn)
        except:
            print(' -- Failed ' + idn)
            pass

