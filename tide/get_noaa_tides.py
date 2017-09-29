"""
Code to TEST getting NOAA tide height data.
"""
import requests
import xml.etree.ElementTree as ET

# for help with ElementTree
# https://docs.python.org/3/library/xml.etree.elementtree.html

a = ('https://tidesandcurrents.noaa.gov/api/datagetter?'
    + 'begin_date=20130101 00:00'
    + '&end_date=20130101 23:00'
    + '&station=9447130'
    + '&product=hourly_height'
    + '&datum=mllw&units=metric&time_zone=gmt'
    + '&application=University of Washington'
    + '&format=xml')
    
b = requests.get(a)
root = ET.fromstring(b.text)

"""
XML Structure
<?xml version="1.0" encoding="UTF-8" ?>
<data>
<metadata id="9447130"  name="Seattle" lat="47.6026" lon="-122.3393"/>
<observations>
<hr t="2013-01-01 00:00"  v="2.467" s="0.020" f="0,0" />
<hr t="2013-01-01 01:00"  v="2.725" s="0.007" f="0,0" />
...
<hr t="2013-01-03 23:00"  v="1.251" s="0.034" f="0,0" />
</observations>
</data>
"""

m_dict = dict()
m = root.find('metadata')
for key in m.keys():
    m_dict[key] = m.attrib[key]

t_list = []
eta_list = []
for e0 in root.findall('observations'):
    for e in e0.findall('hr'):
        t_list.append(e.attrib['t'])
        eta_list.append(float(e.attrib['v']))

