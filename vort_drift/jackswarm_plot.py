"""
First processing and plotting of the Jack Swarm data.
"""
# setup
import os; import sys
alp = os.path.abspath('/Users/pm7/Documents/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
import zfun; reload(zfun)
import matfun; reload(matfun)

import numpy as np
import matplotlib.pyplot as plt

do_movie = True
testing = False

# specification of time loop
dtmin = .5 # time in minutes between movie frames
Dtmin = 5 # time in minutes for data to highlight in each frame
if testing:
    tmin1_arr = np.array([53]) #np.arange(-25+Dtmin,0,dtmin)
else:
    tmin1_arr = np.arange(-25+Dtmin,95,dtmin)

def jack_to_df(jnum):
    # Function to load Jack data.
    dir0 = '/Users/pm7/Documents/tools_data/vortex_drifter/2014_11_JackSwarm_Data/'
    fn = 'JACK' + str(jnum) + '.CSV'
    # load a jack data frame, and make a datetime column
    import pandas as pd
    jdf = pd.read_csv(dir0 + fn, parse_dates = [['Date',' Time']])
    # remove spaces from column headings
    cols = jdf.columns
    cols = cols.map(lambda x: x.replace(' ', '')
        if isinstance(x, (str, unicode)) else x)
    jdf.columns = cols    
    # and specify the index column
    jdf1 = jdf.set_index('Date_Time')        
    return jdf1

jlist = [108, 205, 317, 405, 511]
jind = dict( zip(jlist, range(len(jlist))) )
# NOTE: Jack 601 throws ValueError: 'Date' is not in list,
# and the file format looks different from the others. 

# put all the data frames in a dict
Jdf = dict()
for jnum in jlist:
    Jdf[jnum] = jack_to_df(jnum)

# add columns to each data frame for time (seconds) and vorticity (rad s-1)  
for jnum in jlist:
    jdf = Jdf[jnum]
    dt = jdf.index.values
    # get a single starting time, from the first jack in the list
    if jnum==108:
        dt0 = dt[0]
    
    # make a time axis (seconds from dt0) out of the datetime format index
    ttd = (dt-dt0) # array of timedeltas
    tf = ttd.astype('float') # convert to array of floats
    ts = tf/1e9 # convert to seconds (timedeltas are in nanoseconds (1e-9 sec))
    
    head = jdf['Heading'].values
       
    # vorticity
    ndel = 2 # 2 would be like center differencing
    dts = ts[ndel:] - ts[:-ndel]
    c1 = np.cos(np.deg2rad(head[:-ndel]))
    s1 = np.sin(np.deg2rad(head[:-ndel]))
    c2 = np.cos(np.deg2rad(head[ndel:]))
    s2 = np.sin(np.deg2rad(head[ndel:]))
    dth = np.arcsin(s2*c1 - c2*s1)
    dth = zfun.filt_hanning(dth, n=60)
    vort = - 2 * dth / dts  
    vort = zfun.filt_hanning(vort, n=60)
    vort = np.concatenate((vort, ndel*[np.nan])) # restore length
    # assumes ndel is a small number (like 2) so offset is negligible
    
    # put the new fields into the data frame
    jdf['ts'] = ts
    jdf['vort'] = vort    
    # put it back in the dict of data frames
    Jdf[jnum] = jdf

# prepare to make movie frames    
outdir = '/Users/pm7/Documents/tools_output/pydev_out/js_plots/'
if do_movie:
    Lfun.make_dir(outdir, clean=True)

# initialize holder for the swarm statistics
sw = dict()
sw['ts'] = 60*(tmin1_arr - Dtmin/2)
sw['x'] = np.nan * np.zeros((len(sw['ts']), len(jlist)))
sw['y'] = np.nan * np.zeros((len(sw['ts']), len(jlist)))
sw['u'] = np.nan * np.zeros((len(sw['ts']), len(jlist)))
sw['v'] = np.nan * np.zeros((len(sw['ts']), len(jlist)))
sw['vort'] = np.nan * np.zeros(len(sw['ts']))

pp = 0 # a counter for numbering the output files
for tmin1 in tmin1_arr:
    tmin0 = tmin1 - Dtmin
    # tmin0,1 are the time limits in minutes for the focus view
        
    nouts = ('0000' + str(pp))[-4:]       
    outname = 'js_' + nouts + '.png'
    outfile = outdir + outname

    # PLOTTING
    plt.close()
    fig = plt.figure(figsize=(15,8))
    print 'plotting ' + outname
    
    ms0 = .3 # small marker
    ms1 = 5 # big marker
    gr0 = 3*[.85] # default gray
    
    # MAP
    ax1 = fig.add_subplot(221)
    # coastline
    Ldir = Lfun.Lstart(alp)
    coast_file = Ldir['data'] + 'coast/pnw_coast_combined.mat'
    cmat = matfun.loadmat(coast_file) 
    ax1.plot(cmat['lon'],cmat['lat'], '-k', linewidth=.5)
    # extras
    ll_lims = [-122.91, -122.865, 47.15, 47.165]
    ax1.axis(ll_lims)
    zfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.ticklabel_format(useOffset=False)
    ax1.set_title('Full Deployment')
    # add data, distinguishing mode (marker and size) and Jack # (color)
    modes = range(4) 
    clist = 'rbgmc'
    cdict = dict(zip(jlist, clist))
    for jnum in jlist:
        # Mode is a value between 0 and 3 - (TRANSIT, WAIT, DRIFT, STOP)
        jdf = Jdf[jnum]
        for mode in modes:
            if mode == 0:
                ax1.plot(jdf.ix[jdf['Mode']==mode, 'Longitude'],
                    jdf.ix[jdf['Mode']==mode, 'Latitude'],
                    marker='.', markersize=ms0,
                    color='k', linestyle='None')
            elif mode == 2:
                ax1.plot(jdf.ix[jdf['Mode']==mode, 'Longitude'],
                    jdf.ix[jdf['Mode']==mode, 'Latitude'],
                    marker='.', markersize=ms0,
                    color=cdict[jnum], linestyle='None')
            else:
                ax1.plot(jdf.ix[jdf['Mode']==mode, 'Longitude'],
                    jdf.ix[jdf['Mode']==mode, 'Latitude'],
                    marker='.', markersize=ms0,
                    color=gr0, linestyle='None')
    
    # MOVING FOCUS BOX    
    ax2 = fig.add_subplot(122)        
    # time range
    ts0 = tmin0 * 60
    ts1 = tmin1 * 60
    
    def deg2xy(lon, lat, lon0, lat0):
        # assumes all inputs are degrees
        # output is x, y (meters)
        RE = 6371e3 # Earth radius (m)
        import numpy as np
        lonr = np.deg2rad(lon)
        latr = np.deg2rad(lat)
        lon0r = np.deg2rad(lon0)
        lat0r = np.deg2rad(lat0)
        x = RE*(lonr - lon0r)*np.cos(lat0r)
        y = RE*(latr - lat0r)
        return x, y
        
    def xy2deg(x, y, lon0, lat0):
        # assumes x, y are meters, and lon0, lat0 are degrees
        # output is lon, lat (degrees)
        RE = 6371e3 # Earth radius (m)
        import numpy as np
        lat0r = np.deg2rad(lat0)
        lon = lon0 + np.rad2deg(x/(RE*np.cos(lat0r)))
        lat = lat0 + np.rad2deg(y/RE)
        return lon, lat
    
    # determine center for box (average all Jack locations in focus time period)
    lon0_list = []
    lat0_list = []

    for jnum in jlist:
        jdf = Jdf[jnum]
        jdff = jdf[jdf['ts']>ts0]
        jdff = jdff[jdff['ts'].values<=ts1]
        lon0_list.append(jdff['Longitude'].mean())
        lat0_list.append(jdff['Latitude'].mean())
    lon0 = np.nanmean(lon0_list)    
    lat0 = np.nanmean(lat0_list)
           
    # vorticity scatterplot in the focus box
    Jdff = dict()  
    for jnum in jlist:
        jdf = Jdf[jnum]
        x, y = deg2xy(jdf['Longitude'].values, jdf['Latitude'].values,
            lon0, lat0)
        jdf['x'] = x
        jdf['y'] = y
        
        jdff = jdf[jdf['ts']>ts0]
        jdff = jdff[jdff['ts'].values<=ts1]
        
        Jdff[jnum] = jdff
        
        for mode in modes:
            if mode == 0:
                ax2.plot(jdf.ix[jdf['Mode']==mode, 'x'],
                    jdf.ix[jdf['Mode']==mode, 'y'],
                    marker='.', markersize=ms0,
                    color='k', linestyle='None')
            elif mode == 2:
                ax2.plot(jdf.ix[jdf['Mode']==mode, 'x'],
                    jdf.ix[jdf['Mode']==mode, 'y'],
                    marker='.', markersize=ms0,
                    color=cdict[jnum], linestyle='None')
            else:
                ax2.plot(jdf.ix[jdf['Mode']==mode, 'x'],
                    jdf.ix[jdf['Mode']==mode, 'y'],
                    marker='.', markersize=ms0,
                    color=gr0, linestyle='None')
        

        # plot vorticity in this focus time, only when we are drifting        
        jdff2 =jdff[jdff['Mode']==2]    
        vscale = 70 
        # positive vorticity      
        jdff2p = jdff2[jdff2['vort']>0]
        ax2.scatter(jdff2p['x'],jdff2p['y'],
        s=abs(vscale*jdff2p['vort'].values),
        c='r',
        edgecolor='r',
        alpha=.2)
        # negative vorticity
        jdff2n = jdff2[jdff2['vort']<0]
        ax2.scatter(jdff2n['x'],jdff2n['y'],
        s=abs(vscale*jdff2n['vort'].values),
        c='b',
        edgecolor='b',
        alpha=.2)
        
        # add vorticity scale markers
        ax2.scatter(.05, .95,
            s=vscale*1, c='r', edgecolor='r', alpha=.2,
            transform=ax2.transAxes)
        ax2.text(.07, .95, r'Jack Vorticity $= + 1 s^{-1}$',
            verticalalignment='center',
            transform=ax2.transAxes)
        ax2.scatter(.05, .9,
            s=vscale*1, c='b', edgecolor='b', alpha=.2,
            transform=ax2.transAxes)
        ax2.text(.07, .9, r'Jack Vorticity $= - 1 s^{-1}$',
            verticalalignment='center',
            transform=ax2.transAxes)
               
        # extra settings
        ax2.axis(200*np.array([-1, 1, -1, 1]))
        ax2.set_aspect(1)
        ax2.set_xlabel('X (m)')
        ax2.set_ylabel('Y (m)')
        ax2.set_title('Focus Box')
        
        # we would also like to estimate the vorticity of the swarm.
        # Need to get x, y, u and v over the time of the focus box.       
        # Get velocity by fitting a line to time and position.
        try:
            pf = np.polyfit(jdff['ts'], jdff['x'], 1)
            pf = np.polyfit(jdff['ts'], jdff['y'], 1)
            sw['u'][pp, jind[jnum]] = pf[0]           
            sw['v'][pp, jind[jnum]] = pf[0]            
            sw['x'][pp, jind[jnum]] = np.nanmean(jdff['x'])
            sw['y'][pp, jind[jnum]] = np.nanmean(jdff['y'])
        except:
            sw['u'][pp, jind[jnum]] = np.nan           
            sw['v'][pp, jind[jnum]] = np.nan            
            sw['x'][pp, jind[jnum]] = np.nan
            sw['y'][pp, jind[jnum]] = np.nan
                
   
    # add the focus box to the map plot
    [x1, x2, y1, y2]  = ax2.axis()
    [lon1, lon2], [lat1, lat2] = xy2deg(np.array([x1, x2]), np.array([y1, y2]),
        lon0, lat0)  
    ax1.plot([lon1, lon2, lon2, lon1, lon1], [lat1, lat1, lat2, lat2, lat1],
        '-', color=gr0, linewidth=3)
    
    # TIME SERIES
    ax3 = fig.add_subplot(223)
    for jnum in jlist:
        jdf = Jdf[jnum]
        jdff = Jdff[jnum]
        
        for mode in modes:
            if mode == 0:
                ax3.plot(jdf.ix[jdf['Mode']==mode, 'ts']/60,
                    0*jdf.ix[jdf['Mode']==mode, 'vort'],
                    marker='.', markersize=ms0,
                    color='k', linestyle='None')
            elif mode == 2:
                ax3.plot(jdf.ix[jdf['Mode']==mode, 'ts']/60,
                    jdf.ix[jdf['Mode']==mode, 'vort'],
                    marker='.', markersize=ms0,
                    color=cdict[jnum], linestyle='None')
                ax3.plot(jdff.ix[jdf['Mode']==mode, 'ts']/60,
                    jdff.ix[jdf['Mode']==mode, 'vort'],
                    marker='.', markersize=ms1,
                    color=cdict[jnum], linestyle='None')
            else:
                ax3.plot(jdf.ix[jdf['Mode']==mode, 'ts']/60,
                    0*jdf.ix[jdf['Mode']==mode, 'vort'],
                    marker='.', markersize=ms0,
                    color=gr0, linestyle='None')
                    
                 
        ax3.plot(jdf['ts']/60, jdf['vort'],'-', color=3*[.9], linewidth=.5)
        ax3.plot(jdff['ts']/60, jdff['vort'],'-',color=cdict[jnum])
    ax3.grid()
    ax3.set_xlabel('Time (minutes)')
    ax3.set_ylabel(r'Vorticity $(s^{-1})$')

    # calculate swarm statistics (we are done with the jnum loop) 
    xx = sw['x'][pp, :]
    yy = sw['y'][pp, :]
    uu = sw['u'][pp, :]
    vv = sw['v'][pp, :]
    if not np.isnan(xx).any(): 
        A = np.column_stack((np.ones(len(xx)), xx, yy))
        coeffs_u = np.linalg.lstsq(A, uu)[0]
        coeffs_v = np.linalg.lstsq(A, vv)[0]    
        sw['vort'][pp] = coeffs_v[1] - coeffs_u[2]
    else:
        sw['vort'][pp] = np.nan

    ax3.plot(sw['ts'][:pp+1]/60, 10*sw['vort'][:pp+1], '*k')
    
    ax3.plot(.05, .05, '*k', transform=ax3.transAxes)
    ax3.text(.07, .05, r'$10 \times$ Swarm Vorticity',
        verticalalignment='center',
        transform=ax3.transAxes)

    
           
    if do_movie:        
        plt.savefig(outfile)
    else:
        plt.show()
        
    pp = pp+1

print 'DONE'