# -*- coding: utf-8 -*-
"""
Class for looking at netcdf files.

"""

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

class Nco:
    
    def __init__(self, fn):
        # initialize some fields
        self.fn = fn
        self.ds = nc.Dataset(fn)
        
    def info(self):
        # prints detailed information about all the variables
        for vn in self.ds.variables:            
            print(self.ds.variables[vn])
                
    def vlist(self):
        # prints variable names organized by number of dimensions
        # e.g. for ii == 4, that is 3D+Time varibles like salt
        for ii in range(5):
            print('\n' + 10*'*' + ' ndim = ' + str(ii) + ' ' + 10*'*')
            for vn in self.ds.variables:
                    if len(self.ds.variables[vn].dimensions) == ii:
                        print(vn)
                        
    def plot_slice(self, vn):
        # plot a slice of a variable
    
        if len(self.ds.variables[vn].dimensions) == 4:
            v = self.ds.variables[vn][0, -1, :, :].squeeze()
        elif len(self.ds.variables[vn].dimensions) == 3:
            v = self.ds.variables[vn][0, :, :].squeeze()
        elif len(self.ds.variables[vn].dimensions) == 2:
            v = self.ds.variables[vn][:, :].squeeze()
        else:
            print('Need a 2D or 2D+Time or 3D+Time variable')
            
        jimax = np.unravel_index(np.argmax(v), v.shape)
        vmax = v[jimax[0], jimax[1]]
        jimin = np.unravel_index(np.argmin(v), v.shape)
        vmin = v[jimin[0], jimin[1]]
        
        self.jimax = jimax
        self.jimin = jimin
        
        plt.close()
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cs = ax.pcolormesh(v)
        fig.colorbar(cs)
        ax.plot(jimax[1], jimax[0], '*k')
        ax.plot(jimin[1], jimin[0], 'ok')
        ax.set_title('%s :: max(*) = %0.2f, min(o) = %0.2f' % (vn, vmax, vmin))
        
        plt.show()
    

