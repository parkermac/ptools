function [lon,lat]=wcofs_xy_2_lonlat(x,y,phi0,theta0)

phi=-y;
theta=x;
[lon,lat]=rotgrd_2_lonlat(phi,theta,phi0,theta0);
