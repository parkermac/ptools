This is code for processing and plotting observational weather data.

NCDC is the National Climate Data Center, which has long term weather data.

To get to the daily data go to a place like:

https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND:USW00024233/detail

which is the SeaTac station, with data from 1948-present.

and more generally you can use the search tools here (be sure to select Daily summaries):

https://www.ncdc.noaa.gov/cdo-web/

Some data ranges:
							stations/[ ]/detail
SeaTac:		1948-present	GHCND:USW00024233	1527887.csv (TSUN = 1965-1996)
Seattle:	1891-present	CITY:US530018 (bug in data retrieval)
UW:			1909-1983		USC00457478 (not long enough into present)

Good data to look at: (be sure to ask for custom csv file, metric units, full date range, and then select just these.  You access all these by first adding the station to your Cart, then looking at the item in your cart.)

TMAX Max temp
TMIN Min temp
PRCP Precipitation (0.1 mm per day)
TSUN Total sunshine for the period (minutes)

The wind data doesn't have enough coverage.

The files you get (by email link, after a short wait) are named by you order number, for example the one I just got, for SeaTac, is called 1527887.csv.

================================================================
* plot_weather.py analyzes and plots weather data

Input: ptools_data/ncdc/1527887.csv

Output: Plot of weather variables, focusing on possible long-term changes such as warming, or seasonal timing.  These end up in ptools_output/ncdc.