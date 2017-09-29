"""
Code, mainly from Doug Latournell, for getting Canadian Tide Data.

"""
import requests

def get_all_perm_dfo_wlev(start_date, end_date):
    """Get water level data for all permanent DFO water level sites
    for specified period.

    :arg start_date: Start date; e.g. '01-JAN-2010'.
    :type start_date: str

    :arg end_date: End date; e.g. '31-JAN-2010'
    :type end_date: str

    :returns: Saves text files with water level data at each site
    """
    stations = {
        'Point Atkinson': 7795,
        'Vancouver': 7735,
        'Patricia Bay': 7277,
        'Victoria Harbour': 7120,
        'Bamfield': 8545,
        'Tofino': 8615,
        'Winter Harbour': 8735,
        'Port Hardy': 8408,
        'Campbell River': 8074,
        'New Westminster': 7654,
    }
    for ttt in stations:
        get_dfo_wlev(stations[ttt], start_date, end_date)


def get_dfo_wlev(station_no, start_date, end_date, outdir, year):
    """Download water level data from DFO site for one DFO station
    for specified period.

    :arg station_no: Station number e.g. 7795.
    :type station_no: int

    :arg start_date: Start date; e.g. '01-JAN-2010'.
    :type start_date: str

    :arg end_date: End date; e.g. '31-JAN-2010'
    :type end_date: str

    :returns: Saves text file in outdir with water level data at one station
    """
    # Name the output file
    outfile = outdir + 'tide_'+str(station_no)+'_'+str(year)+'.csv'
    
    # Form urls and html information
    base_url = 'http://www.meds-sdmm.dfo-mpo.gc.ca/isdm-gdsi/twl-mne/inventory-inventaire/'
    form_handler = (
        'data-donnees-eng.asp?user=isdm-gdsi&region=PAC&tst=1&no='
        + str(station_no))
    sitedata = {
        'start_period': start_date,
        'end_period': end_date,
        'resolution': 'h',
        'time_zone': 'u', # 'l' for local, 'u' for UTC
    }
    data_provider = (
        'download-telecharger.asp'
        '?File=E:%5Ciusr_tmpfiles%5CTWL%5C'
        + str(station_no) + '-'+start_date + '_slev.csv'
        '&Name=' + str(station_no) + '-'+start_date+'_slev.csv')
    # Go get the data from the DFO site
    with requests.Session() as s:
        s.post(base_url + form_handler, data=sitedata)
        r = s.get(base_url + data_provider)
    # Write the data to a text file
    with open(outfile, 'w') as f:
        f.write(r.text)
        

station_dict = {
    'Point Atkinson': 7795,
    'Vancouver': 7735,
    'Patricia Bay': 7277,
    'Victoria Harbour': 7120,
    'Bamfield': 8545,
    'Tofino': 8615,
    'Winter Harbour': 8735,
    'Port Hardy': 8408,
    'Campbell River': 8074,
    'New Westminster': 7654,
    }
    
# the program to run
import os
outdir = os.environ.get('HOME') + '/Documents/ptools_data/tide/'

# note that the start and end dates are Local (Pacific) time, so we need to adjust to get
# a clean year in UTC, and we can only go by days
year = 2013
for sn in station_dict.keys():
    if sn == 'Point Atkinson':
        print('Working on ' + sn)
        get_dfo_wlev(station_dict[sn], '31-DEC-' + str(year-1), '01-JAN-' + str(year+1), outdir, year)

