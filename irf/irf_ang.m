function [d]=irf_ang(t1,p1,t2,p2)
%IRF_ANG   Compute the angle between two vectors in degrees
%
% [d]=irf_ang(t1,p1,t2,p2)
% [d]=irf_ang(v1,v2) % column vectors in radial coordinates
% [d]=irf_ang(v1,v2,whatever) % v1,v2 column vectors in xyz
% returns the angle between two vectors in degrees
% v - column vectors
% t,p given in degrees
% theta corresponds to latitude
% pi is positive against clock, zero at  x-axis
%

if nargin ==3 % two vectors in XYZ
  nt=size(t1,2);
  if (nt>3)
    time=t1(:,1);
  else
    time=[];
  end
  d=[time acos(irf_dot(irf_norm(t1),irf_norm(p1),1))*180/pi];
else
  if nargin ==2  % angles supplied as two radial vectors
    t1r=pi/180*t1(2);
    p1r=pi/180*t1(3);
    t2r=pi/180*p1(2);
    p2r=pi/180*p1(3);
  else
    t1r=pi/180*t1;
    p1r=pi/180*p1;
    t2r=pi/180*t2;
    p2r=pi/180*p2;
  end
  n1=[cos(t1r).*cos(p1r) cos(t1r).*sin(p1r) sin(t1r)];
  n2=[cos(t2r).*cos(p2r) cos(t2r).*sin(p2r) sin(t2r)];
  
  d=acos(dot(n1,n2))*180/pi;
end

