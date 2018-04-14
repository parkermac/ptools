"""
Code to analyze the aestus output in APE(rho,z) space.

"""
import os
import sys
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
import zrfun
import Lfun
Ldir = Lfun.Lstart('aestus1', 'A1')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

dir0 = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

fn_neap = dir0 + 'f2013.03.12/ocean_his_0001.nc'
fn_spring = dir0 + 'f2013.03.19/ocean_his_0001.nc'

G, S, T = zrfun.get_basic_info(fn_neap)
lon_vec = G['lon_rho'][0,:].squeeze()
lat_vec = G['lat_rho'][:,0].squeeze()
ii0, ii1, ifr = zfun.get_interpolant(np.array([0.02, 1.5]), lon_vec)
jj0, jj1, jfr = zfun.get_interpolant(np.array([44.9, 45.1]), lat_vec)
i0 = ii0[0]
i1 = ii1[1]
j0 = jj0[0]
j1 = jj1[1]
h = G['h'][j0:j1, i0:i1]
dx = G['DX'][j0:j1, i0:i1]
dy = G['DY'][j0:j1, i0:i1]
mask = dy = G['mask_rho'][j0:j1, i0:i1]
da = dx*dy
da = np.ma.masked_where(mask==0, da)

ny, nx = da.shape
zr, zw = zrfun.get_z(h, 0*h, S)
dzr = np.diff(zw, axis=0)
dv = dzr * da.reshape(1,ny,nx)

# APE contours
z_edges = np.linspace(-20,0,21)
rho_edges = np.linspace(0,27,28)
z_bins = z_edges[:-1] + np.diff(z_edges)/2
rho_bins = rho_edges[:-1] + np.diff(rho_edges)/2
Ze, Re = np.meshgrid(z_edges, rho_edges)
Zb, Rb = np.meshgrid(z_bins, rho_bins)
g = 9.8
APE = g * (Re - 27) * Ze

# plotting
plt.close('all')
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,12))

pp = 0
for fn in [fn_neap, fn_spring]:

    ds = nc.Dataset(fn)
    
    # APE
    rho = ds['rho'][0, :, j0:j1, i0:i1].squeeze()
    ape = g * (rho - 27) * zr
    ape_ave = np.sum(ape*dv) / np.sum(dv)
    rho_ave = np.sum(rho*dv) / np.sum(dv)
    
    print('rho_ave = %0.2f' % (rho_ave))
    
    # salinity variance
    s = ds['salt'][0, :, j0:j1, i0:i1].squeeze()
    sbar = np.sum(s*dv) / np.sum(dv)
    svar = (s - sbar)**2
    
    print('ape_ave = %0.2f J m-3' % (ape_ave))
    ds.close()
    
    ax = axes[0,pp]    
    cs = ax.contour(Re, Ze, APE, np.linspace(0, 6000, 13))
    plt.clabel(cs, inline=1, fontsize=10, fmt='%d')
    
    for ii in range(nx):
        for jj in range(ny):
            ax.plot(rho[:,jj,ii], zr[:,jj,ii], alpha=.2)    
    ax.set_xlim(0,27)
    ax.set_ylim(-20,0)
    ax.set_xlabel('Sigma (kg m-3)')
    ax.grid(True)
    if pp==0:
        ax.set_title('Neap (red line = avg. Sigma, contours = APE)')
        ax.set_ylabel('Z (m)')
    elif pp==1:
        ax.set_title('Spring')
    # add a line at the average density
    ax.plot([rho_ave,rho_ave], [-20,0], '-r', linewidth=3)
        
    ax = axes[1,pp]
    
    svar_fresh = svar[s<=sbar]
    svar_salty = svar[s>sbar]
    ape_fresh = ape[s<=sbar]
    ape_salty = ape[s>sbar]
    
    ax.plot(svar_fresh.flatten(), ape_fresh.flatten(), '.b',
            markersize=.3, label='fresher than average')
    ax.plot(svar_salty.flatten(), ape_salty.flatten(), '.r',
            markersize=.3, label='saltier than average')
    #ax.plot(svar.flatten(), ape.flatten(), '.r', markersize=1)
    
    ax.grid(True)
    ax.set_xlim(0,500)
    ax.set_ylim(0,2100)
    ax.set_xlabel('Salinity Variance (g/kg)^2')
    if pp==0:
        ax.set_ylabel('APE (J m-3)')
        ax.legend()
    
    
    
    pp += 1

plt.show()