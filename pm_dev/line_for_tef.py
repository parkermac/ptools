"""
Code to develop an algorithm for making diagonal lines for TEF extractions
and river carving.

It creates the indices for a stair-step path between two points (which are
also defined by indices on the grid).  It always takes a step of one grid index
either in the x or y direction, never diagonally.

For a TEF section these would be interpreted as indices on the psi-grid, and then
you would have to pull out u- or v-grid transports depending on the direction of
the step.

If this is for carving a river channel these would be interpreted as indices on the
rho-grid.

RESULT: it works perfectly.

"""

import matplotlib.pyplot as plt
import numpy as np

# make a grid of indices (integers)
N = 20
x = np.arange(N)
y = np.arange(N)
X,Y = np.meshgrid(x,y)

# initialize plotting
plt.close('all')
fig = plt.figure(figsize=(14,10))

def dist(x, x1, y, y1):
    # standard distance between two points
    d = np.sqrt((x1-x)**2 + (y1-y)**2)
    return d
    
def dist_normal(x0, x1, y0, y1, xp, yp):
    # This finds the shortest distance from xp, yp to the segment
    # defined by points (x0, y0) and (x1, y1).
    
    # make sure x is increasing (inside the function)
    if x1 < x0:
        (x0, x1) = (x1, x0)
        (y0, y1) = (y1, y0)
    dx = x1-x0
    dy = y1-y0
    
    # find xp2, yp2: the position of the point on the line
    # that is closest to xp, yp
    if dy == 0: # line is horizontal
        xp2 = xp
        yp2 = y0
    elif dx == 0: # line is vertical
        yp2 = yp
        xp2 = x0
    else: # The line is sloping: then we find the equation for
        # a line through xp, yp that is normal to the original line,
        # and then find the point xp2, yp2 where the two lines intersect.
        m = dy/dx
        b = y0 - m*x0
        bp = yp + (1/m)*xp
        xp2 = (bp - b)/(m + (1/m))
        yp2 = m*xp2 + b
    # and return the distance between the point xp, yp and the intersection point
    dxp = xp2 - xp
    dyp = yp2 - yp
    dist_n = np.sqrt(dxp**2 + dyp**2)
    return dist_n
    
def get_stairstep(x0, x1, y0, y1):
    # brute force method of generating stairstep path: always choosing
    # the point with the shortest distance to the end will generate
    # the straightest path.
    d = dist(x0, x1, y0, y1)
    xx = []
    yy = []
    xx.append(x0)
    yy.append(y0)
    x = x0
    y = y0
    while d > 0:
                
        # find distances of all 4 surrounding points that
        # are allowed choices to the end of line
        dn = dist(x, x1, y+1, y1)
        ds = dist(x, x1, y-1, y1)
        de = dist(x+1, x1, y, y1)
        dw = dist(x-1, x1, y, y1)
        # close_arr is an array that is positive for surrounding points
        # that are closer to the endpoint than the current point
        close_arr = np.array([d-dn, d-ds, d-de, d-dw])
        
        # find distances of all 4 surrounding points that
        # are allowed choices to closest point on line
        dn_n = dist_normal(x0, x1, y0, y1, x, y+1)
        ds_n = dist_normal(x0, x1, y0, y1, x, y-1)
        de_n = dist_normal(x0, x1, y0, y1, x+1, y)
        dw_n = dist_normal(x0, x1, y0, y1, x-1, y)
        # gather them in an array
        norm_arr = np.array([dn_n, ds_n, de_n, dw_n])
        
        # step_array is an array of the steps we take to make each
        # of the surrounding points
        step_arr = np.array([[0,1],[0,-1],[1,0],[-1,0]])
        
        # closer is a Boolean array that is True for all step choices
        # that got us closer to the target
        closer = close_arr >= 0
        
        # make shorter versions of step_arr and norm_arr that only include
        # choices that got us closer
        step_arr = step_arr[closer]
        norm_arr = norm_arr[closer]
        
        # then find the index of the remaining choice that was closest to the line
        imin = np.argmin(norm_arr)
        
        # take the step in the chosen direction
        step = step_arr[imin,:]
        x += step[0]
        y += step[1]
        xx.append(x)
        yy.append(y)
        
        # and update the distance so the loop knows when it is done
        d = dist(x, x1, y, y1)
        
    # pack results as arrays
    XX = np.array(xx)
    YY = np.array(yy)
    
    # check that result never steps more or less than one
    DD = np.abs(np.diff(XX)) + np.abs(np.diff(YY))
    if (DD==1).all():
        pass # the result passes this test
    else:
        print('Error in result stepsize')
    
    # check that endpoints are correct
    if (XX[0] != x0) or (XX[-1] != x1) or (YY[0] != y0) or (YY[-1] != y1):
        print('Error in result endpoints!')
        
    return XX, YY
    
# loop over a number of different starting and ending points,
# which are indices on the grid
for case in range(6):
    Nn = N-4
    # Vary Nn from 2 to N-1 to change the amount
    # of variation of the index that varies less.
    # The other index always spans the full range: 0 to N-1.
    if case == 0:
        # flatter than 1:1
        x0, x1 = (0, N-1)
        y0, y1 = (2, Nn)
    elif case == 1:
        # flatter than 1:1 negative
        x0, x1 = (0, N-1)
        y0, y1 = (Nn, 2)
    elif case == 2:
        # flatter than 1:1 negative, reverse order
        x0, x1 = (N-1, 0)
        y0, y1 = (2, Nn)
    elif case ==3:
        # steeper then 1:1
        x0, x1 = (2, Nn)
        y0, y1 = (0, N-1)
    elif case ==4:
        # steeper then 1:1 negative
        x0, x1 = (2, Nn)
        y0, y1 = (N-1, 0)
    elif case ==5:
        # steeper then 1:1 negative, reverse order
        x0, x1 = (Nn, 2)
        y0, y1 = (0, N-1)
    # get the solution
    XX, YY = get_stairstep(x0, x1, y0, y1)
    # PLOTTING
    ax = fig.add_subplot(2,3,case+1)
    # draw the grid
    for ii in range(N):
        ax.plot(X[ii,:],Y[ii,:],'-k', alpha=.5)
        ax.plot(X[:,ii],Y[:,ii],'-k', alpha=.5)
    # defining line
    ax.plot([x0, x1], [y0, y1], '*-y', ms=24, lw=3)
    ax.plot(x0, y0, '*r', ms=24) # starting point
    # the stairstep solution
    ax.plot(XX,YY, '-og', lw=3)
plt.show()
