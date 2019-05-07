clear all;

% Cutout indices for the subsample:
ii=237:390; 
jj=714:886;

HomeDir='/mnt/lfs3/projects/nosofs/Alexander.Kurapov/WCOFS/';


grdfile=[HomeDir 'Prm/grd_wcofs_large_visc200.nc'];
N=40;
Vtransform=2;
Vstretching=4;
theta_s=8;     % parameter for stretching near surface
theta_b=3;     % parameter for stretching near bottom
Tcline=50;     % thermocline depth

h=ncread(grdfile,'h');

h=h(ii,jj);

[nx,ny]=size(h);

% s-coord of the section:
zeta0=zeros(nx,ny);
kgrid=0; % 0 for rho, 1 for W points
column=1;
plt=0;
z_r=zeros(nx,ny,N);
for jsect=1:ny
 [z,sc,Cs]=scoord_new(h,zeta0,theta_s,theta_b,Tcline,N,...
                     kgrid,column,jsect,Vtransform,Vstretching,plt);
 z_r(:,jsect,:)=reshape(z,[nx 1 N]);
end

z_u=0.5*(z_r(1:end-1,:,:)+z_r(2:end,:,:));
z_v=0.5*(z_r(:,1:end-1,:)+z_r(:,2:end,:));
