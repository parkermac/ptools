This is for processing the Ecology CTD and Bottle data.

======================================================================
* process_station_info.py makes a DateFrame of station locations, etc.

Input: ptools_data/ecology/raw/ParkerMacCreadyCoreStationInfoFeb2018.xlsx

Output: ptools_data/ecology/sta_df.p with data like
         Desig                                            Descrip     ...       Latitude   Longitude
Station                                                               ...                           
ADM001       C                         Admiralty Inlet - Bush Pt.     ...      48.029813 -122.617933

and the full column listing is
['Desig', 'Descrip', 'Basin', 'Max_Depth', 'Latitude', 'Longitude']

and the Basins have values:
{'Admiralty Inlet',
 'Grays Harbor',
 'Hood Canal Basin',
 'Main Basin',
 'South Basin',
 'Strait of Georgia',
 'Strait of Juan de Fuca',
 'Whidbey Basin',
 'Willapa Bay'}
 
======================================================================
* station_list.py prints a screen listing of the stations for use in javascript for a webpage map.  Entries are like:
sta_list = [
    {sta:"ADM001", lat:48.0298, lng:-122.6179},
	
======================================================================
* station_map.py plots station locations from sta_df.p on two maps, one for the Salish Sea and one for the coastal estuaries.

======================================================================
* process_casts.py does initial processing of the CTD cast data

Input: ptools_data/ecology/raw/ParkerMacCready2017CTDDataFeb2018.xlsx and
	   ptools_data/ecology/raw/ParkerMacCready1999-2016CTDDataMay2018.xlsx
	which I had gotten by sepcial request to Ecology (Julia Bos, or Skip Albertson)

Output: ptools_data/ecology/Casts_[1999-2017].p with data like:
        Salinity  Temperature      Sigma       Chl        DO      Turb      Z Station       Date
0      30.129299       8.1719  23.427200  0.426924  8.512160  0.789888 -125.5  ADM001 2012-01-11
1      30.122299       8.1754  23.421200  0.432654  8.492985  0.739318 -125.0  ADM001 2012-01-11
2      30.113501       8.1798  23.413700  0.452711  8.488438  0.728413 -124.5  ADM001 2012-01-11

Issue: what are the units of DO and chl?  Answer (from plot_ctd.py):
data_names =  ['Salinity','Temperature','Sigma', 'Chl', 'DO',   'Turb', 'Z']
data_units =  ['psu',     'deg C',      'kg/m3', 'ug/l', 'mg/l', '',    'm']

======================================================================
* process_bottles.py does initial processing of the Bottle data

Input: ptools_data/ecology/raw/Parker_2006-present_Nutrients.xlsx (goes 2006-2017)

Output: ptools_data/ecology/Bottles_[2006-2017].p with data like:
     Station       Date       Z  PO4(uM)D  SiOH4(uM)D   NO3(uM)D  NO2(uM)D  NH4(uM)D        DIN  Znom
0     ADM001 2006-02-06 -27.323  2.305970   55.271759  26.999115  0.129438  0.072036  27.200589 -30.0
1     ADM001 2006-02-06 -10.131  2.337134   56.988705  27.365324  0.124822  0.055482  27.545627 -10.0
2     ADM001 2006-02-06  -1.077  2.378515   58.525520  27.768095  0.120208  0.083661  27.971964  -0.0
3     ADM001 2006-03-21 -30.813  2.238526   53.678001  26.457401  0.272849  0.455226  27.185476 -30.0
4     ADM001 2006-03-21 -10.786  2.171502   50.757725  25.796255  0.253642  0.518785  26.568682 -10.0
5     ADM001 2006-03-21  -0.931  2.047256   48.457275  24.136850  0.229608  0.633126  24.999585  -0.0

======================================================================
* plot_ctd_casts.py makes nice plots CTD casts (property vs. z, colored by month) at each station.  Because this can also be used with LiveOcean model output, it is not modified to be part of the SSMSP workflow.

Input: ptools_data/ecology/Casts_2017.p & ptools_data/canada/Casts_2017.p
	   and model extractions like
	   LiveOcean_output/cast/cas4_v2_lo6biom/PSB003_2017.12.11.nc
	   
Output: plots, one per station, names like ADM001.png, in folders like
	ptools_output/ecology/casts/ (used in a lab exercise for Coastal 320 class in 2018), or
	ptools_output/ecology/casts_cas4_v2_lo6biom/ if you include the model
	
======================================================================
* process_obsmod_series.py makes DataFrames comparing obs & mod at 3-4 depths.  Because this can also be used with LiveOcean model output, it is not modified to be part of the SSMSP workflow.

Input: ptools_data/ecology/Casts_2017.p & ptools_data/canada/Casts_2017.p &
	   ptools_data/ecology/Bottles_2017.p
	   and model extractions like
	   LiveOcean_output/cast/cas6_v3_lo8b/PSB003_2017.12.11.nc
	   [or other years like 2018 and 2019]
	   
Output: ptools_output/ecology/ObsMod_cas6_v6_lo8b_2017.p (used in the series and scatter plots below) with data like
	
     Station       Date  Znom         ...          Mod DIN (uM) Density (kg m-3) Mod Density (kg m-3)
0     ADM001 2017-02-23     0         ...               18.3425          1022.36              1021.87
1     ADM001 2017-02-23   -10         ...               20.3044          1022.59              1022.91
2     ADM001 2017-02-23   -30         ...               20.8001          1022.67              1023.41
3     ADM001 2017-04-28     0         ...               8.61608          1019.62              1020.71
4     ADM001 2017-04-28   -10         ...               14.4416           1021.7              1022.27
5     ADM001 2017-04-28   -30         ...               17.2305          1022.16              1022.84

the full list of columns is:
['Station', 'Date', 'Znom', 'Z', 'Salinity', 'Temp. (deg C)',
       'Chl (mg m-3)', 'DO (mg L-1)', 'DIN (uM)', 'Mod Salinity',
       'Mod Temp. (deg C)', 'Mod Chl (mg m-3)', 'Mod DO (mg L-1)',
       'Mod DIN (uM)', 'Density (kg m-3)', 'Mod Density (kg m-3)']
	   
======================================================================
* plot_obsmod_series.py makes time series plots comparing obs & mod over a year (2017,8,9).  Because this can also be used with LiveOcean model output, it is not modified to be part of the SSMSP workflow.

Input: ptools_output/ecology/ObsMod_cas6_v6_lo8b_2017.p
	   
Output: plots, one per station, names like ADM001.png, in folders like
	ptools_output/ecology/val_series_cas4_v2_lo6biom/ or
	ptools_output/ecology/web_series_cas4_v2_lo6biom/ (formatted for the website, e.g. no map)
	

======================================================================
* plot_obsmod_scatter.py makes time a scatterplot comparing obs & mod over a year.  Because this can also be used with LiveOcean model output, it is not modified to be part of the SSMSP workflow.

Input: ptools_output/ecology/ObsMod_cas6_v3_lo8b_2017.p
	   
Output: a scatterplot of all data (for three depths)

======================================================================
* plot_obsmod_map.py [OBSOLETE] makes plots that are map summary info about fields at three depths, comparing modeled and observed fields.  This makes circles to compare obs-mod, with area proportional to value.  I decided I did not like them, and went to version _2 below.

Input: ptools_output/ecology/ObsMod_cas6_v3_lo8b_2017.p

Output: ptools_output/ecology/obsmod_maps/cas6_v3_lo8b_2017_DIN.png for example

======================================================================
* plot_obsmod_map_2.py makes plots that are map summary info about fields at three depths, comparing modeled and observed fields.  This one uses colored side-by-side boxes for the comparison.

Input: ptools_output/ecology/ObsMod_cas6_v3_lo8b_2017.p

Output: ptools_output/ecology/obsmod_maps/cas6_v3_lo8b_2017_DIN.png for example

======================================================================
* bin_by_region.py replicates a step in obs_collias.  It gathers all the casts in a Region, for all years, and then for each cast interpolates in the vertical to just have data at [-30, -10, 0] m.

Input: ptools_data/ecology/Casts_[1999-2017].p &
	   ptools_data/ecology/Bottles_[2006-2017].p

Output: ptools_output/ecology/region_[1-8].py which are DataFrames containing data like:
           Station  Z (m) Salinity Temp. (deg C) DO (mg L-1) NO3 (uM)
Date                                                                 
2009-03-12  RSR837      0  30.3119        7.0192     8.41159  27.2884
2009-03-12  RSR837    -10  30.3119         7.017     8.41581      NaN
2009-03-12  RSR837    -30  30.3116        7.0192     8.42661      NaN
2009-04-07  RSR837      0  30.3256        7.5435     8.56143  23.6588
2009-04-07  RSR837    -10  30.3283        7.4309     8.60177      NaN
2009-04-07  RSR837    -30  30.3287        7.4321     8.65763      NaN

