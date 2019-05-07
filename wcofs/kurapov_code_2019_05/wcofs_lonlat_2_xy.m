function [x,y]=wcofs_lonlat_2_xy(lon,lat,phi0,theta0,regridYES)

% Regular 2D arrays of "local" coordinates x and y of WCOFS, which is 
% a spherical grid in rotated coordinates with the pole at phi0,theta0
% x increases along the local latitude
% y increases in the direction opposite to local longitude
%
% if regridYES==1 then
% x and y are regularized using meshgrid to eliminate roundup errors and
% make x and y suitable to interp2. 

[phi,theta]=lonlat_2_rotgrd(lon,lat,phi0,theta0);
if regridYES
 [y,x]=meshgrid(-phi(1,:),theta(:,1));
else
 y=-phi;
 x=theta;
end
