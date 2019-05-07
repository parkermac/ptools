function [phi,theta]=rotgrd_2_lonlat(phi2,theta2,phi0,theta0)

phi2=phi2*pi/180;
theta2=theta2*pi/180;
phi0=phi0*pi/180;
theta0=theta0*pi/180;

phi=phi0+atan2(sin(phi2).*cos(theta2),...
           cos(phi2).*cos(theta2).*sin(theta0)+sin(theta2).*cos(theta0));
theta=asin(-cos(phi2).*cos(theta2).*cos(theta0)+sin(theta2).*sin(theta0));

phi=phi*180/pi;
theta=theta*180/pi;
