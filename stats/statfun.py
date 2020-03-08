# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 09:55:33 2016

@author: PM5

Functions for statistics.

"""

import numpy as np
from scipy import stats

def cross_correlation(x,y):
    """
    Returns the cross-correlation coefficient for a pair of time series, x
    and y.
    It is defined as in Oke et al. Part I (2002) JPO, vol 32, p. 1379,
    Eqn. (C1)
    """
    # first flatten the records
    x = x.flatten()
    y = y.flatten()
    # then remove all nans
    mask = (x!=np.nan) & (y!=np.nan)
    x = x[mask]
    y = y[mask]
    if len(x)>0 and len(y)>0:
        cc = ( (x-x.mean()) * (y-y.mean()) ).mean() / (x.std() * y.std()) 
    else:
        cc = np.nan
    return cc
    
def trend_test(x,y):
    '''
    Test code for fitting a line to some points, including the calculation
    of statistics and their confidence intervals
    
    Note: stats.t.ppf(1-.025,100) = matlab tinv(1-.025,100)    
    '''
    # make some data to fit
    N = 20
    x = np.linspace(0,1,N)
    y_0 = 0
    dydx = 1
    y_std = .1
    y_noise = y_std * np.random.randn(1,N)
    y_nonoise = y_0 + dydx * x
    y = y_0 + dydx * x + y_noise

    % fit the data
    Y = y.T
    X = np.ones(N,1), x.reshape((N,1))]
    B = np.inv(X.T.dot(X)).dot( (X.T.dot(Y)) )
    yy_0 = B(1)
    dyydx = B(2)
    yy = yy_0 + dyydx * x

    # % fit the data using "polyfit"
    # % NOTE: p(1) = dydx, p(2) = y(x=0);
    # p = polyfit(x,y,1);
    # yy_fit = polyval(p,x);

    # calculate the 95% confidence interval for the mean and trend
    # NOTE the call to "tinv" for the inverse cumulative Student's t-distribution
    # s_eps = sqrt( sum((y-yy).^2)/(N-2) );
    # s_x = std(x);
    # ci_mean = std(y) * tinv(1-0.05/2,N-2) / sqrt(N);
    # ci_trend = s_eps * tinv(1-0.05/2,N-2) / sqrt((N-1)*s_x*s_x);
    #
    # figure
    # Z_fig(14)
    # plot(x,y_nonoise,'-b',x,y,'*b',x,yy,'-or',x,yy_fit,'-g')
    # legend('Data - no noise','Data','Fit','Polyfit',0);
    # [xt,yt] = Z_lab('ul');
    # text(xt,yt,['dydx = ',num2str(dydx)],'color','b');
    # [xt,yt] = Z_lab('ll');
    # text(xt,yt,['dydx = ',num2str(dyydx),' \pm ',num2str(ci_trend)],'color','r');
    # [xt,yt] = Z_lab('ur');
    # text(xt,yt,['mean = ',num2str(mean(y)),' \pm ',num2str(ci_mean)],'color','r', ...
    #     'horizontalalignment','r');
    # title(['real std = ',num2str(y_std),', s_{\epsilon} = ',num2str(s_eps)]);
    
