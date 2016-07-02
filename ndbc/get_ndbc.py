# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 08:06:36 2016

@author: PM5

Code to automate getting multiple years of NDBC buoy data
"""

import bs4
import urllib.request as U

# Specify which NDBC buoy number
#sn = 46022 # near Eel River
#sn = 46029 # near Coumbia River
#sn = 42040 # near Mississippi River

sn = 'wpow1' # Near West Point

yr_list = range(1984,2016)

for yr in yr_list:
    try:
        idn = str(sn) + 'h' + str(yr)
        print('Attempting to get ' + idn)
        url_str = ('http://www.ndbc.noaa.gov/view_text_file.php?filename=' +
                   idn + '.txt.gz&dir=data/historical/stdmet/')
        html = U.urlopen(url_str, timeout=10)
        soup = bs4.BeautifulSoup(html, 'html.parser')
        sn_text = soup.findAll(text=True)
        sns = str(sn_text)[2:-2]
        sns = sns.replace('\\n','\n')
        f = open(idn + '.txt','w')
        f.write(sns)
        f.close()
    except:
        pass


