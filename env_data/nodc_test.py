"""
Code to get data from an NODC buoy.
"""

# import some modules 
import urllib2
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

# define the station
sta = 46029
#sta = 46050

# define timestrings for the time limits of the request
t1 = datetime.now()
ndays = 5
t0 = t1 - timedelta(ndays) # ndays ago (is 5 the most we can request?)
t0_str = datetime.strftime(t0, '%Y-%m-%d') + 'T00:00Z'
t1_str = datetime.strftime(t1, '%Y-%m-%d') + 'T00:00Z'

# createt the url string
url_str = ('http://sdf.ndbc.noaa.gov/sos/server.php'
    + '?request=GetObservation&service=SOS&version=1.0.0'
    + '&offering=urn:ioos:station:wmo:' + str(sta)
    + '&observedproperty=sea_water_temperature'
    + '&responseformat=text/xml;schema=%22ioos/0.6.1%22'
    + '&eventtime=' + t0_str +'/' + t1_str)

# get the XML tree and root
try:
    file = urllib2.urlopen(url_str, timeout = 10)
    tree = ET.parse(file)
    root = tree.getroot()
except:
    print 'problem downloading XML'
 
# initialize empty lists to store the data
data = []
taxis = []

# get the data                                                                          
try:
    rt = root.tag
    for e0 in root.findall(".//"):
        if 'Quantity' in e0.tag:
            if e0.attrib['name'] == 'WaterTemperature':
                data.append(float(e0.text))
                data_units = e0.attrib['uom']
        elif 'timePosition' in e0.tag:
            taxis.append(datetime.strptime(e0.text, '%Y-%m-%dT%H:%M:%SZ'))
        elif 'StationName' in e0.tag:
            sname = e0.text
except:
    print 'problem parsing data from XML'

# PLOTTING

import matplotlib.pyplot as plt
plt.close()

# create an hour axis, just for plotting
hours = []
for ta in taxis:
    tdel = ta - taxis[0]
    hours.append(tdel.total_seconds()/ 3600.)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.plot(hours, data, '-og')
ax.set_xlabel('Hours from '
    + datetime.strftime(taxis[0], '%Y-%m-%dT%H:%M:%SZ'))
ax.set_ylabel('Temperature deg ' + data_units)
ax.set_title('Station #' + str(sta) + ': ' + sname)

plt.show()