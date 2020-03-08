import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# make some data to fit
N = 1000
x = np.linspace(0,1,N)
y_0 = 0
dydx = 1
y_std = .1
y_noise = y_std * np.random.randn(N)
y = y_0 + dydx * x + y_noise

# fit the data
Y = y.reshape(N,1)
X = np.concatenate( (np.ones((N,1)), x.reshape(N,1)), axis=1)
B = np.linalg.inv(X.T.dot(X)).dot( (X.T.dot(Y)) )
yy_0 = B[0]
dyydx = B[1]
yy = yy_0 + dyydx * x


# calculate the 95% confidence interval for the mean and trend
# NOTE the call to "tinv" for the inverse cumulative Student's t-distribution
s_eps = np.sqrt( ((y-yy)**2/(N-2)).sum() );
s_x = np.std(x);
# Note: stats.t.ppf(1-.025,100) = matlab tinv(1-.025,100)
ci_mean = np.std(y) * stats.t.ppf(1-0.05/2,N-2) / np.sqrt(N);
ci_trend = s_eps * stats.t.ppf(1-0.05/2,N-2) / np.sqrt((N-1)*s_x*s_x);

fs = 14

#plt.close('all')
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ax.plot(x, y, '*b', label='Data')
ax.plot(x, yy, '-r', label='Fit Line')
lh = ax.legend(loc='lower right')
#plt.setp(lh.get_texts(), fontsize=fs) 
plt.setp(ax.get_legend().get_texts(), fontsize=fs) 

ax.text(.05, .9, '$Mean = %0.2f \pm %0.3f$' % (y.mean(), ci_mean),
    transform=ax.transAxes, size=fs)
ax.text(.05, .8, '$Slope = %0.2f \pm %0.3f$' % (dyydx, ci_trend),
    transform=ax.transAxes, size=fs)
ax.text(.05, .7, '95% Confidence Intervals',
    transform=ax.transAxes, size=fs)
    
ax.set_title('N = %d' % (N), size=fs+2)
    
plt.show()
