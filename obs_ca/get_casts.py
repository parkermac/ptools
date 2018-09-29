"""
Code to automate getting Canadian CTD and Bottle data.

NOTE: this code does not currently work.  I need to figure out how
to get into the password protected site.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import requests

# first load a spreadsheet of cast info
indir = Ldir['parent'] + 'ptools_data/canada/'
fn = indir + 'IOS_Search_Water_Profiles_Data_27_42_2017.csv'
df = pd.read_csv(fn, header=0)
castname = []
for item in df.FILENAME:
    castname.append(item.split('/')[-1])    
df['castname'] = castname
df = df.set_index('castname')

print('Code not working yet')

#cfn = df.index[0]
# #This URL will be the URL that your login form points to with the "action" tag.
# POST_LOGIN_URL = 'https://www.waterproperties.ca'
# #This URL is the page you actually want to pull down with requests.
# REQUEST_URL = 'https://www.waterproperties.ca/osd_data_archive/' + cfn
# payload = {
#     'username': 'parker',
#     'pass': 'pogoCaRiv'
# }
# with requests.Session() as session:
#     post = session.post(POST_LOGIN_URL, auth=('parker','pogoCaRiv'))
#     r = session.get(REQUEST_URL)
#     print(r.text)   #or whatever else you want to do with the request data!
 
        
# for cfn in [df.index[0]]:
#     print('-- getting: ' + cfn)
#     url = 'https://www.waterproperties.ca/osd_data_archive/' + cfn
#     response = requests.get(url, auth=('parker','pogoCaRiv'))
#     outfile = indir + 'raw/' + cfn + '.txt'
#     with open(outfile, mode='wb') as localfile:
#         localfile.write(response.content)
#     response.close()


