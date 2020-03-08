% makes a 3d plot of cool LO stuff

clear
close all

phome = '/Users/pm7/Documents/';
indir = 'LiveOcean_roms/output/cas6_v3_lo8b/f2019.07.04/';
infile = [phome,indir,'ocean_his_0001.nc'];
lon = nc_varget(infile,'lon_rho');
lat = nc_varget(infile,'lat_rho');

eta = squeeze(nc_varget(infile,'zeta'));
z = -squeeze(nc_varget(infile,'h'));
eta = squeeze(nc_varget(infile,'zeta'));
DO = squeeze(nc_varget(infile,'oxygen',[0, 0, 0, 0],[1, 1, -1, -1]));

%
figure;
set(gcf,'position',[100 100 1000.5 600]);

%h1 = surfl(lon,lat,eta);

h2 = surfl(lon,lat,z);

aa = axis;

axis([-130 -122 42 52 -3000 3])
shading flat
colormap copper
lighting phong
%set(gca,'DataAspectRatio',[1.46628 1 5]);
set(gca,'view',[55 34]);
set(gca,'fontsize',14);

%light('Position',[-125, -65, 10000],'Style','local')

