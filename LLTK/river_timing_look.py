"""
Code to look at changes in TIMING of Skagit River flow.

"""

import pandas as pd
import matplotlib.pyplot as plt
in_dir = '/Users/PM5/Documents/LiveOcean_data/rivers/Data_historical/'

plt.close('all')
fig = plt.figure(figsize=(12,8))

rn_list = ['fraser', 'skagit', 'nisqually', 'duckabush']

counter = 1
for rn in rn_list:
    
    ax = fig.add_subplot(2,2,counter)
    
    in_file = in_dir + rn + '.p'

    q = pd.read_pickle(in_file)

    qq = q.resample('m', loffset='-15d').mean()

    yrs = set(qq.index.year)

    df = pd.DataFrame(index=range(1,13), columns=yrs)

    for yr in yrs:
        qy = qq[qq.index.year == yr]
        qy.index = qy.index.month
        qs = pd.Series(index=range(1,13))
        # we introduce qs to allow for missing values
        qs[qy.index] = qy.values    
        df[yr] = qs.values
        # df: index = month, columns = year

    dft = df.transpose()
    # dft: index = year, columns = month

    sdft = dft[[4,5,6,7,8]] # just use selected months
    
    sdftc = sdft.cumsum(axis=1) * (30*86400/1e9)

    # plot cumulative sum versus year
    use_legend = False
    if counter==4:
        use_legend = True
    sdftc.plot(ax=ax, linewidth=3, legend=use_legend, xlim=(1980,2015), 
        marker='o', title='Cumulative flow (km3): ' + rn.title())
        
    counter += 1

plt.show()
