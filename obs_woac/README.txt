Programs for processing WOAC cruise data.  These are CTD+Bottle casts of data provided by Simone Alin (11/2018).  The code is structured to be much like that in obs_ecy and obs_collias.

======================================================================
* process_data.py

Input: ptools_data/woac/raw/WOAC_data_9-7-2018_dataToParker.xlsx

Output: ptools_data/woac/sta_df.p
		ptools_data/woac/Casts_2017.p
		ptools_output/woac/station_map.png
		
sta_df is a DataFrame with info like:
        Latitude Longitude
Station                   
W5       47.8799  -122.372
W7       47.9831  -122.622
W28      47.7054  -122.453

Casts is a DataFrame with data like:
      Cruise       Date       Time  Longitude  ...    Omega Ar     DO (uM)       Z (m)  ncast
0    CAB1065 2017-04-06   21:14:35 -122.37086  ...    0.703018  228.296664 -240.105631      0
1    CAB1065 2017-04-06   21:18:06 -122.37211  ...         NaN  229.989032 -150.736001      0
2    CAB1065 2017-04-06   21:19:43 -122.37268  ...    0.801193  254.515428 -110.611965      0
3    CAB1065 2017-04-06   21:21:21 -122.37311  ...         NaN  264.230826  -81.203427      0
4    CAB1065 2017-04-06   21:22:59 -122.37269  ...         NaN  264.173253  -50.635715      0
5    CAB1065 2017-04-06   21:24:08 -122.37247  ...         NaN  272.537131  -30.371035      0
6    CAB1065 2017-04-06   21:25:06 -122.37276  ...         NaN  279.368144  -20.883421      0
7    CAB1065 2017-04-06   21:24:59 -122.37271  ...    0.971966  279.066720  -20.763447      0
8    CAB1065 2017-04-06   21:26:19 -122.37327  ...    1.433686  332.598505  -11.275388      0
9    CAB1065 2017-04-06   21:27:20 -122.37362  ...         NaN  386.030887   -5.545058      0
10   CAB1065 2017-04-06   21:28:15 -122.37372  ...         NaN  396.535737   -2.485974      0
11   CAB1065 2017-04-06   21:28:08 -122.37374  ...    1.943827  396.939987   -2.461184      0
12   CAB1065 2017-04-06   23:11:33 -122.62151  ...    1.068248  270.851023  -90.494133      1
13   CAB1065 2017-04-06   23:13:58 -122.62246  ...         NaN  275.819606  -50.747273      1
14   CAB1065 2017-04-06   23:15:45 -122.62286  ...         NaN  285.464407  -30.306311      1
and the full list of columns is:
['Cruise', 'Date', 'Time', 'Longitude', 'Latitude', 'Station',
       'Pressure (dbar)', 'Temp. (deg C)', 'Salinity', 'Sigma (kg m-3)',
       'DO Flag', 'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'Chl (ug L-1)',
       'TA1 (umol kg-1)', 'DIC1 (umol kg-1)', 'TA1 Flag', 'DIC1 Flag',
       'TA2 (umol kg-1)', 'DIC2 (umol kg-1)', 'TA2 Flag', 'DIC2 Flag', 'pH',
       'pCO2 (uatm)', 'CO2 (umol kg-1)', 'HCO3- (umol kg-1)',
       'CO3-- (umol kg-1)', 'Omega Ca', 'Omega Ar', 'DO (uM)', 'Z (m)', 'ncast']

Issue: because I define casts as having distinct Stations and Dates, it is possible that some casts will be cut into two parts becasue they extended over UTC midnight.  Does this matter?