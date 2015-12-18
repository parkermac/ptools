"""
Test of zfun.get_interpolant()

Parker MacCready
"""

# setup
import os; import sys
alp = os.path.abspath('../../LiveOcean/alpha/')
if alp not in sys.path:
    sys.path.append(alp)
import zfun; reload(zfun)

import numpy as np

# INPUT - both must be numpy arrays
#
# x is the array of points we want interpolants for
x = np.array([2.2, 3.5])
# x = 2.2 # test of non-array input
#
# xvec is the coordinate axis (must be monotonic and increasing)
xvec = np.arange(10)
# xvec[5] = -3 # test of non-monotonic xvec

# run function
b = zfun.get_interpolant(x, xvec)

for bb in b:
    print ''
    print 'ind0 = ' + str(bb[0])
    print 'ind1 = ' + str(bb[1])
    print 'a    = ' + str(bb[2])
