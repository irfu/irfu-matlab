function [yo,zo,tilt]=irf_gse2gsm(yr,day,sec,yi,zi,ig);
%IRF_GSE2GSM change from GSE to GSM
%
%   [yo,zo,tilt]=irf_gse2gsm(yr,day,sec,yi,zi,ig);
%   change from GSE to GSM coordinate
%   yr is a 2-digit year, day is day no (Jan 1=Day 1)
%   sec is sec of day (it can be an array same dim as yse zse)
%   yi,zi are the input y and z coordinates
%   ig indicates coord conversion direction (0=GSE to GSM)
%   yo,zo are the output y and z coordinates
%   tilt is the dipole tilt angle: +ve towards sun
%
% $Id$

%	calculate constants and dipole tilt angle
fdoy=sec/86400;
d=day+365.25*yr-0.75+fdoy;
vl0=4.881628+0.0172027912*d;
g=6.256584+0.0172019697*d;
ecc=0.0167504484-1.1444e-9*d;
vl=vl0-9.92394e-5+2*ecc.*sin(g).*(1+1.25*ecc.*cos(g));
tau=vl0+2*pi*fdoy;
sl=sin(vl); cl=cos(vl);
st=sin(tau); ct=cos(tau);
x1=-.1860153*st-.0685841*ct;
x2=.0685841*st-.1860153*ct;
smu=cl.*x1+sl.*(.3978291*.9801505-.9174595*x2);
tilt=180*asin(smu)/pi;
d2=cl.*(.9801505*.3978291-.9174595*x2)-sl.*x1;
d3=.9801505*.9174595+.3978291*x2;
dp=sqrt(d2.*d2+d3.*d3);
d2=d2./dp; d3=d3./dp;

%	do the coordinate transformation
if ig==0	% convert from GSE to GSM
	yo=d3.*yi-d2.*zi;
	zo=d2.*yi+d3.*zi;
else
	yo=d3.*yi+d2.*zi;
	zo=-d2.*yi+d3.*zi;
end

