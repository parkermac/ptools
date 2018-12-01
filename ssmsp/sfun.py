"""
This module contains utility functions for interpolation, filtering
and inspection.

Parker MacCready
"""

import netCDF4 as nc
import numpy as np

import os 
from datetime import datetime, timedelta
import csv
import shutil

def make_dir(dirname, clean=False):
    # Make a directory if it does not exist.
    # Use clean=True to clobber the existing directory.
    if clean == True:
        import shutil
        shutil.rmtree(dirname, ignore_errors=True)
        os.mkdir(dirname)
    else:
        try:
            os.mkdir(dirname)
        except OSError:
            pass # assume OSError was raised because directory already exists
            
# ========================================================================

def get_interpolant(x, xvec, extrap_nan=False):
    """
    Returns info to allow fast interpolation.

    Input:
    x = data position(s) [1-D numpy array]
    xvec = coordinate vector [1-D numpy array without nans]
        NOTE: xvec must be monotonically increasing

    *kwargs*
    Set extrap_nan=True to return nan for the fraction
    whenever an x value is outside of the range of xvec.
    The default, with extrap_nan=False, is that
    if the x is out of the range of xvec, or if x=nan it returns
    the interpolant for the first or last point.
    E.g. [0, 1, 0.] for x < xvec.min()
    
    Output: three 1-D numpy arrays of the same size as x
    i0 = index below [int]
    i1 = index above [int]
    fr = fraction [float]

    If the x is ON a point in xvec the default is to return
    the index of that point and the one above, with fr=0,
    unless it is the last point in which case it is the index
    of that point and the point below, with fr = 1.
    """
    
    from warnings import filterwarnings
    filterwarnings('ignore') # skip some warning messages
    
    def itp_err(message='hi'):
        print('WARNING from get_interpolant(): ' + message)

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

    return i0, i1, fr

def find_nearest(array, value):
    # gives the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_ind(array, value):
    # gives the index of the item in array that is closest to value
    idx = (np.abs(array-value)).argmin()
    return idx

def filt_hanning(data, n=40):
    """
    Input: 1D numpy array data
    Output: Array of the same size, filtered with Hanning window of length n,
        padded with nan's
    If n=1 it just returns matrix you gave it
    """
    if n == 1:
        smooth = data
    else:
        filt = hanning_shape(n=n)
        npad = np.floor(len(filt)/2).astype(int)
        smooth = np.convolve(data, filt, mode = 'same')
        smooth[:npad] = np.nan
        smooth[-npad:] = np.nan
    return smooth
    
def filt_hanning_mat(data, n=40):
    """
    Input: ND numpy array, with time on axis 0.
    Output: Array of the same size, filtered with Hanning window of length n,
        padded with nan's
    """
    filt = hanning_shape(n=n)
    filt = filt / filt.sum()
    n = np.ceil(len(filt)/2).astype(int)
    sh = data.shape
    df = data.flatten('F')
    dfs = np.convolve(df, filt, mode = 'same')
    smooth = dfs.reshape(sh, order='F')
    smooth[:n,:] = np.nan
    smooth[-n:,:] = np.nan
    return smooth

def filt_godin(data):
    """
    Input: 1D numpy array of HOURLY values
    Output: Array of the same size, filtered with 24-24-25 Godin filter,
        padded with nan's
    """
    filt = godin_shape()
    filt = filt / filt.sum()
    n = np.ceil(len(filt)/2).astype(int)
    smooth = np.convolve(data, filt, mode = 'same')
    smooth[:n] = np.nan
    smooth[-n:] = np.nan
    return smooth

def filt_godin_mat(data):
    """
    Input: ND numpy array of HOURLY, with time on axis 0.
    Output: Array of the same size, filtered with 24-24-25 Godin filter,
        padded with nan's
    """
    filt = godin_shape()
    filt = filt / filt.sum()
    n = np.ceil(len(filt)/2).astype(int)
    sh = data.shape
    df = data.flatten('F')
    dfs = np.convolve(df, filt, mode = 'same')
    smooth = dfs.reshape(sh, order='F')
    smooth[:n,:] = np.nan
    smooth[-n:,:] = np.nan
    return smooth
    
def godin_shape():
    """
    Based on matlab code of 4/8/2013  Parker MacCready
    Returns a 71 element numpy array that is the weights
    for the Godin 24-24-25 tildal averaging filter. This is the shape given in
    Emery and Thomson (1997) Eqn. (5.10.37)
    ** use ONLY with hourly data! **
    """
    k = np.arange(12)
    filt = np.NaN * np.ones(71)
    filt[35:47] = (0.5/(24*24*25))*(1200-(12-k)*(13-k)-(12+k)*(13+k))
    k = np.arange(12,36)
    filt[47:71] = (0.5/(24*24*25))*(36-k)*(37-k)
    filt[:35] = filt[:35:-1]
    return filt

def hanning_shape(n=40):
    """
    Returns a Hanning window of the specified length.
    """
    ff = np.cos(np.linspace(-np.pi,np.pi,n+2))[1:-1]
    filt = (1 + ff)/2
    filt = filt / filt.sum()
    return filt

def earth_rad(lat_deg):
    """
    Calculate the Earth radius (m) at a latitude
    (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid

    INPUT: latitude in degrees

    OUTPUT: Earth radius (m) at that latitute
    """
    a = 6378.137 * 1000; # equatorial radius (m)
    b = 6356.7523 * 1000; # polar radius (m)
    cl = np.cos(np.pi*lat_deg/180)
    sl = np.sin(np.pi*lat_deg/180)
    RE = np.sqrt(((a*a*cl)**2 + (b*b*sl)**2) / ((a*cl)**2 + (b*sl)**2))
    return RE

def ll2xy(lon, lat, lon0, lat0):
    # This converts lon, lat into meters relative to lon0, lat0.
    # It should work for lon, lat scalars or arrays.
    # NOTE: lat and lon are in degrees!!
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    x = R * clat * np.pi * (lon - lon0) / 180
    y = R * np.pi * (lat - lat0) / 180
    return x, y
    
def get_rc(NP):
    # figure out near-optimal numer of rows and columns for plotting
    NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
    NC = np.ceil(NP/NR).astype(int)
    return NR, NC
    
def get_irc(ii, NC):
    # get row and column of plot ii when there are NC columns
    ir = int(np.floor(ii/NC))
    ic = int(ii - NC*ir)
    return ir, ic






