function [ed,d]=irf_edb(e,b,angle_lim)
%IRF_EDB   Computer Ez under assumption E.B=0
%
% [ed,d]=irf_edb(e,b,angle_lim)
% ed - E field with assumption E.B=0
%  d - B elevation angle above spin plane
%  
%  e - [Ex Ey] or [t Ex Ey] or [t Ex Ey Ez] (Ez whatever)
%  b - [Bx By Bz] or [t Bx By Bz]
%  angle_lim - B angle with respec to the spin plane should be at least angle_lim degrees
%             otherwise Ez is put to 0
%
% $Id$

global AV_DEBUG 
if isempty(AV_DEBUG), debug=0; else debug=AV_DEBUG;end

if nargin==0, help irf_edb;return;end
if nargin == 2, 
  angle_lim=20; 
  if debug == 1, disp('Using limiting angle of 20 degrees');end
end

le= size(e,2); % 
lb= size(b,2); % 

if size(b,1) ~= size(e,1),
 if (debug == 1), disp('E and B are not of the same length. Interpolating B.');end
 b=irf_resamp(b,e);
end


if le < 2
 disp('E has not enough components');return;
elseif le == 2
 ed=[e(:,1) e(:,2) e(:,1)*0];
elseif le == 3
 ed=[e(:,2) e(:,3) e(:,1)*0];
else
 ed=[e(:,2) e(:,3) e(:,1)*0];
end

if lb == 2
 disp('B has not enough components');return;
elseif lb == 3
 bd=[b(:,1) b(:,2) b(:,3)];
elseif lb > 3
 bd=[b(:,2) b(:,3) b(:,4)];
end


d=atan2(bd(:,3),sqrt(bd(:,1).^2+bd(:,2).^2))*180/pi;
ind=find(abs(d)>angle_lim);

if ~(size(ind)==[0 0])
ed(ind,3)=-(ed(ind,1).*bd(ind,1)+ed(ind,2).*bd(ind,2))./bd(ind,3);
end

tt=sprintf('E.B=0 was used for %0.0f from %0.0f data points',length(ind),length(d));
if debug == 1, disp(tt);end

if le > 2
 ed=[e(:,1) ed];
end
