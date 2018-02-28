#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 09:44:40 2018

@author: pm7

Code to test scraping data from the Ecology Marine Water Quality Site.

The bit I want to scrape for latitude and longitude:
    
<tr><td  class="data" style="padding-top:8;">
<span class="datalabel">station information:</span><br>

	<table width="396" style="text-align:center;font-size:70%" cellpadding="3" cellspacing="0" border="1">
	<tr style="color:#666666"><td>type</td><td>latitude</td><td>longitude</td><td>max depth</td><td>county</td><td>wria</td><td>*map detail <!-- <span style="color:#cc6600">(external&nbsp;site)</span> --></td></tr>
	<tr style="background-color:#fdfdfd;text-align:center"><td>core</td><td>48.0300</td><td>-122.6167</td><td>153 meters</td><td>Island</td><td>06</td><td><a class="googlelink" target="_blank" href="http://maps.google.com/maps?q=48.03+-122.6167+(ADM001)">Google maps&reg;</a>&nbsp;</td></tr>
	</table>
</td></tr>


"""

import requests
import bs4
import pandas as pd

# where the data csv files are, and where to save the DataFrame
dir0 = '/Users/pm7/Documents/ptools_data/ecology/'
fn_out = dir0 + 'sta_df.p'

# set up a DataFrame to save the results in
sta_df = pd.DataFrame(columns=['ID','longitude','latitude'])

for staID_try in range(1,10):#200):
    
    url_str = 'https://fortress.wa.gov/ecy/eap/marinewq/mwdataset.asp?staID='+str(staID_try)
    
    DATA_URL = 'http://wateroffice.ec.gc.ca/report/real_time_e.html'
    response = requests.get(url_str)
    #print(response.status_code)
    soup = bs4.BeautifulSoup(response.content, 'lxml')    
    data = soup.find_all('td', class_='data')
    
    for item in data:        
        ii = item.find_all(name='input')
        for iii in ii:    
            try:
                if iii.attrs['name'] == 'staname':
                    staname = iii.attrs['value']
            except KeyError:
                pass           
            try:
                if iii.attrs['name'] == 'staID':
                    staID = iii.attrs['value']
            except KeyError:
                pass
        
        if 'station information' in item.get_text():           
            sinfo = item           
            rows = sinfo.find_all('tr')            
            nn = 0
            for col in rows[0].find_all('td'):
                if col.get_text() == 'latitude':
                    nlat = nn
                    #print(nlat)
                elif col.get_text() == 'longitude':
                    nlon = nn
                    #print(nlon)           
                nn += 1    
            nn = 0
            for col in rows[1].find_all('td'):
                if nn == nlat:
                    latitude = float(col.get_text())
                    #print(col.get_text())
                elif nn == nlon:
                    longitude = float(col.get_text())
                    #print(col.get_text())           
                nn += 1 
    if str(staID_try) == staID and len(staname) > 0:
        
        sta_df.loc[staname,['ID','longitude','latitude']] = [staID, longitude, latitude]
        
        print('%s %s %0.3f %0.3f' % (staname, staID, longitude, latitude)) 

sta_df.to_pickle(fn_out)            
