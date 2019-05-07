% AK, 2019/5/01: Demo to sample a profile from a cutout made for Parker:

clear all;
dir0 = '/Users/pm7/Documents/LiveOcean_roms/output/WCOFS_avg_Exp37/';
fname= [dir0, 'zuvts_ORWA_Parker_Exp37_0010.nc'];
grdfile= [dir0, 'grd_wcofs_large_visc200.nc'];

% Mooring coordinates (NH10):
lato=44.65;
lono=-124.3;

ii=375:696;
jj=1030:1710;

% The north pole of the rotated spherical system:
theta0=37.4;
phi0=-57.6;

N=40;

% END OF INPUTS ^^^^^^^^^^^^^^^^
nx=length(ii);
ny=length(jj);


lon_rho=ncread(grdfile,'lon_rho');
lat_rho=ncread(grdfile,'lat_rho');
mask_rho=ncread(grdfile,'mask_rho');
h=ncread(grdfile,'h');

lon_rho=lon_rho(ii,jj);
lat_rho=lat_rho(ii,jj);
mask_rho=mask_rho(ii,jj);
h=h(ii,jj);

figure;
hold on;
[cc,hh]=contour(lon_rho,lat_rho,h,[200 1000],'k-');
[cc,hh]=contour(lon_rho,lat_rho,mask_rho,[0.5 0.5],'k-');
set(hh,'linewidth',2);
plot(lono,lato,'ro');

% note: last input variable is 1, meaning that x and y will be regridded 
% to enable interpolation:
[x_rho,y_rho]=wcofs_lonlat_2_xy(lon_rho,lat_rho,phi0,theta0,1);

% note: last input variable is 0 for the obs locations, 
% the first two arguments can be either scalar, or vector, or array
[xo,yo]=wcofs_lonlat_2_xy(lono,lato,phi0,theta0,0);

% Read T from the cutout file:
T=ncread(fname,'temp');
T=double(T);

To=zeros(N,1);
for k=1:N
 To(k)=interp2(y_rho,x_rho,T(:,:,k),yo,xo);
end

figure;
plot(To,[1:N]');






