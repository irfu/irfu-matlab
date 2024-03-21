function [rtp]=car2sph(xyz,direction_flag)
%%
% [rtp] = IRF.CAR2SPH(xyz,direction_flag)
% returns magnitude, theta and phi angle from column vector xyz (first coloumn is x ....)
% theta is 0 at equator
% direction_flag = -1  -> to make transformation in opposite direction


time_flag=0;
if length(xyz(1,:))>3
  time_flag=1;
  time = xyz(:,1);
  xyz=xyz(:,[2 3 4]);
end

if nargin>1
  if direction_flag==-1
    r=xyz(:,1);st=sin(xyz(:,2)*pi/180);ct=cos(xyz(:,2)*pi/180);sp=sin(xyz(:,3)*pi/180);cp=cos(xyz(:,3)*pi/180);
    z=r.*st;
    x=r.*ct.*cp;
    y=r.*ct.*sp;

    if time_flag == 1
      rtp=[time x y z];
    else
      rtp=[x y z];
    end
    return
  end
end

xy=xyz(:,1).^2+xyz(:,2).^2;
r=sqrt(xy+xyz(:,3).^2);
t=atan2(xyz(:,3),sqrt(xy))*180/pi;
p=atan2(xyz(:,2),xyz(:,1))*180/pi;

if time_flag == 1
  rtp=[time r t p];
else
  rtp=[r t p];
end

if (size(r,1) < 2) && nargout == 0
  sprintf('r,theta,phi=%g, %g, %g',r,t,p);
end
