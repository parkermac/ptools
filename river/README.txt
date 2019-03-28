README for river code

This code gets river data from USGS and Environment Canada sources.

===================================================================
* make_river_info.py

This takes various files of river information in LiveOcean_data/rivers/Info_curated/
and creates the master list:
- ptools_output/river/pnw_all_2016_07/river_info.csv
which has entries like:

rname,usgs,ec,nws,ratio,depth,width,max_dist
clowhom,,08GB013,,32.3873,10.0,2.0,0.2
squamish,,08GA022,,2.3449,5.0,2.0,0.2
fraser,,08MF005,,1.0848,10.0,5.0,1.0

for 47 rivers.

It also makes the tracks files:
- ptools_output/river/pnw_all_2016_07/tracks/[river name].csv
with entries like:

ind,lon,lat
0,-124.0919946022139,44.43017917325204
1,-124.0849466259081,44.4232896732934
2,-124.0805232621173,44.42122487204401
3,-124.0758039139337,44.42250877878311

which go from offshore to onshore.

===================================================================
* make_historical.py

This goes through the list of rivers and their gauges and scaling factors here:
ptools_output/river/pnw_all_2016_07/river_info.csv
and creates long data records for each, saved as pickle files in:
LiveOcean_data/rivers/Data_historical/[river name].p
These are pandas Series with entries like:
1980-01-01 12:00:00    124.589141
1980-01-02 12:00:00    103.424960
1980-01-03 12:00:00     93.042532
1980-01-04 12:00:00     84.257400
1980-01-05 12:00:00    148.548592
(units are m3/s and scaling factors have been applied)

These are used by e.g. LiveOcean/forcing/riv2/make_forcing_main.py to fill in
when the requested time period spans the saved historical record.  Otherwise
we go straight to USGS, NWS, or EC.

Currently these go 1980.01.01 through 2018.12.31 (code run 2019.03.20), and
the log and a plot from this run are included in the data folder.

It makes use of the same code which LiveOcean uses to make rivers for
hindcasts and forecasts: LiveOcean/forcing/riv2/river_class.py.

Note: To get Canadian rivers for 2017 and 2018 it uses values scraped from
a LiveOcean run which were made at a time when all were available:
LiveOcean_output/river/cas4_v2_2017.01.01_2018.12.31.p.
We do this because historical values are only available on the web as
tables that we scrape up through the end of 2016.  The "current" values
are only available back 18 months or so, so there is a 6 month gap.
Apparently when I made the forcing files for cas4_v2 it covered this gap,
and then new ones were filled in when forecasts started around October 2018.
There are notes on this in Evernote - LiveOcean - Rivers New.

Note: This version of 2019.03.21 is the first in which I have combined data
from the N. and S. forks of the Skokomish to get a better total Skokomish
(which otherwise cut off most storm flows after 2010).  A similar fix
applies to the Hamma Hamma (see notes below).  This fix is hard-coded
into the method get_usgs_data() in LiveOcean/forcing/riv2/river_class.py.

Note: The skokomish and hamma fixes work for both historical data and for when
we get usgs data for run forcing.  But when we are doing a forecast the
Skokomish has an NWS number and so it does not use the fix.  However, it appears
that what comes from NWS does not have the problem that the main Skok
gauge has, so this should be fine. [need to check to be sure]

Note: The Hamma Hamma was originally created as 0.4459
times the Skokomish flow, which has the problem we just solved above.  To fix
this I used something similar to what I did for the Skokomish.
Note that the scale factor for the original Skok is 1.081, and 0.4459/1.081 = 0.4125,
so I think 0.4125 is the scale factor I should use after constructing the Skokomish
using the 2 gauge method, because this gives the value presumably AFTER multiplying
by 1.081, whereas the 0.4459 factor for the Hamma was to be applied to the Skok
maine gauge BEFORE multiplying by 1.081.

NOTE: There is also a non-negligible flow from the Cushman Dam (peaks around 70 m3/s)
that I should include as a new point source.  The data and notes are in
ptools_data/Cushman_Dam.

Log: I re-ran on 2019.03.21 just for hamma and skokomish.

===================================================================
* test_ec.py

Code to test the timespan covered by EC river data sources, and attempts
at error handling.  One issue is that I found this error in river_class.get_ec_data():
# NEW to deal with the problem that when you request date from before
# 18 months ago it gives the most recent data instead.
and in test_ec.py I implemented a flag to (i) warn that this was happening, and
(ii) to not use that data.

On 2019.03.20 I implemented this fix in riv2/river_class.py.  Since the current cas4
forecast uses riv2 I will need to keep an eye on things when I push the riv2 code to boiler.

Currently the problem should not come up because I now have historical EC data through the end
of 2018, and if I make riv2 forcing files the live EC data they get will be from just the
first few months of 2019.  But by late 2020 if I am trying to make riv2 files for early 2019
it will be a problem.

===================================================================
* make_clim.py

Goes through all the historical files created by make_historical.py and
creates climatologies with a yearday time axis, in:
LiveOcean_data/rivers/Data_clim/[rivername].csv
with entries like:
1,132.02877242289205
2,127.48091381153615
3,114.55168454959363
4,116.94430187268989
5,112.3731493513734
...
363,155.44137040091366
364,164.31967779732648
365,159.53777085255692
366,176.83405361520875
so that we cover leap year.

And I re-ran this on 2019.03.21 to match the work on Data_historical above.

===================================================================
* make_T_historical.py

This gets temperature records where they are available, from 1980 onwards
for USGS, and for the last 18 months for EC rivers.  In construction it is much
like make_historical.py, but the data end up in:
LiveOcean_data/rivers/Data_T_historical/[rivername].p
with are pandas Series with entries like:
1980-01-01 12:00:00     8.000000
1980-01-02 12:00:00     7.750000
1980-01-03 12:00:00     7.750000
1980-01-04 12:00:00     7.750000
(units are Degrees C)

It only has files for 13 rivers.  Currenly it refers to functions in
LiveOcean/forcing/riv1 instead of riv2, but they should be identical.
Also the time periods hard-wired into the code were current a few years
ago and would need to be edited to run again, but for the time being
I see no need to re-run.

===================================================================
* make_T_clim.py

This is analogous to make_clim.py above, and makes .csv files in
LiveOcean_data/rivers/Data_T_clim.  These are the data that are
always used by LiveOcean/forcing/riv2/make_forcing_main.py, e.g.,
to decide temperature.  Often this involves using a proxy river
that has temperature data.

===================================================================
* make_river_infoA.py

This makes the files:
- ptools_output/river/analytical/river_info.csv
- ptools_output/river/analytical/tracks/creek1.csv

* make_QT_climA.py

These make the files:
- LiveOcean_data/rivers/Data_clim/creek1.csv
- LiveOcean_data/rivers/Data_T_clim/creek1.csv

which are analogous to those we use for realistic runs, but for idealized
runs like aestus1.
===================================================================



