function [phi2,theta2]=lonlat_2_rotgrd(phi,theta,phi0,theta0)

phi=phi*pi/180;
theta=theta*pi/180;
phi0=phi0*pi/180;
theta0=theta0*pi/180;

phi1=phi-phi0;
phi2=atan2(sin(phi1).*cos(theta),...
           cos(phi1).*cos(theta).*sin(theta0)-sin(theta).*cos(theta0));
theta2=asin(cos(phi1).*cos(theta).*cos(theta0)+sin(theta).*sin(theta0));

phi2=phi2*180/pi;
theta2=theta2*180/pi;
