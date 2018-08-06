"""
Code to test zfun.get_interpolant()

RESULT: in order to get the expected behaviour when x is equal to
the last point in xvec, AND extrap_nan = True, we added these lines near the end:

            # override for the case where x = the last point of xvec
            fr[X[:,0]==XVEC[0,-1]] = 1.0

and this change is now implemented in zfun.  8/1/2018

"""

import numpy as np
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
def itp_err(message='hi'):
    print('WARNING from get_interpolant(): ' + message)

# inputs
xvec = np.arange(10)

# scan through all important cases
for extrap_nan in [False, True]:
    print(' ')
    for xx in [-1, 0, 5.5, 9, 11]:
        
        x = np.array([xx])

        # input error checking
        if isinstance(x, np.ndarray) and isinstance(xvec, np.ndarray):
            pass # input type is good
        else:
            itp_err('Inputs must be numpy arrays')

        # some preconditioning of the input
        x = x.flatten()
        xvec = xvec.flatten()

        # more error checking
        if np.isnan(x).any() and extrap_nan==False:
            itp_err('nan found in x')
        if np.isnan(xvec).any():
            itp_err('nan found in xvec')
        if not np.all(np.diff(xvec) > 0):
            itp_err('xvec must be monotonic and increasing')

        nx = len(x)
        nxvec = len(xvec)

        X = x.reshape(nx, 1) # column vector
        xvec = xvec.reshape(1, nxvec)
        XVEC = xvec.repeat(nx, axis=0) # matrix

        # preallocate results arrays
        i0 = np.zeros(nx, dtype=int)
        i1 = np.zeros(nx, dtype=int)
        fr = np.zeros(nx, dtype=float)

        # calculate index columns
        mask = X >= XVEC
        # the above line broadcasts correctly even if nx = nxvec
        # because we forced X to be a column vector
        i0 = mask.sum(axis=1) - 1

        # these masks are used to handle values of x beyond the range of xvec
        lomask = i0 < 0
        himask = i0 > nxvec - 2
        i0[lomask] = 0
        i0[himask] = nxvec - 2
        i1 = i0 + 1

        # compute the fraction
        xvec0 = xvec[0,i0]
        xvec1 = xvec[0,i1]
        fr = (x - xvec0)/(xvec1 - xvec0)

        # fractions for out of range x
        if extrap_nan == False:
            fr[lomask] = 0.
            fr[himask] = 1.
        elif extrap_nan == True:
            fr[lomask] = np.nan
            fr[himask] = np.nan
            # override for the case where x = the last point of xvec
            fr[X[:,0]==XVEC[0,-1]] = 1.0
        
        print('x=%0.1f: i0=%d, i1=%d, fr=%0.1f, extrap_nan=%s' % (x, i0, i1, fr, str(extrap_nan)))
