function [E]=irf_e_vxb(v,b,flag)
%IRF_E_VXB   Compute VxB and ExB/B^2
%
% [E]=irf_e_vxb(v,b)
% [V]=irf_e_vxb(e,b,-1)
%
% calculates electric field or velocity E=-vxB or v=ExB/B^2
% v[km/s], B[nT], E[mV/m]
% v=[t vx vy vz] (t..vz column vectors)
% b=[t bx by bz]
% v and b can be at different sampling
% E is at b sampling
% V is at e sampling
%
% $Id$

global AV_DEBUG;if isempty(AV_DEBUG), debug=0; else debug=AV_DEBUG;end

if (nargin ==3) & (flag == -1),
 e=v;
 if size(b,1) ~= size(e,1), 
  if debug == 1, disp('interpolating b to e');end
  bb=c_resamp(b,e);b=bb;clear bb; 
 end
 v=irf_vec_x_scal(irf_cross(e,b),[b(:,1) irf_abs(b,1)],-2);
 v=irf_tappl(v,'*1e3');
 E=v;
else,
  if size(v,1) == 1,    
    if size(v,2)==3
      v=[b(1,1) v];
    elseif size(v,2)<3,
      error('v has to few components');
    end
    v=irf_resamp(v,b); 
  else,                 
    b=irf_resamp(b,v);
  end
  % estimating E =(v x B) 
  E=irf_tappl(irf_cross(v,b),'*1e3*1e-9*1e3*(-1)');
end

