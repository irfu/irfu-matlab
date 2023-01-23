function [emp,nl,nm,nn]=irf_eb_nrf(e,b,v,flag)
%IRF_EB_NRF Find E and B in MP system given B and MP normal vector
%
%       [xmp]=irf_eb_nrf(x,b,v,flag)
% [xmp,L,M,N]=irf_eb_nrf(x,b,v,flag)
%
% xmp=[time x_L x_M x_N]
%
% A) find x vector (ex E) in MP system given B and MP normal vector
%
%  L - along B
%  N - closest to v
%  M - NxL
% the coordinate system follows B and thus is not stationary
%
% B) If flag==1 then find xmp in stationary reference frame defined
%  N - along v
%  L - the mean direction of B in plane perpendicular to N
%  M - NxL
%
% C) If flag==L_vector then find xmp in stationary reference frame defined
%  N - along v
%  L - closest to the direction specified by L_vector (ex: maximum variance direction)
%  M - NxL
%
%  x b  - 4 columns, first time
%  v = [vx vy vz]
%  xmp=[t xl xm xn]
%

if nargin ==3, flag_case='A';end
if (nargin ==4)
  if length(flag)==3
    L_direction=flag;clear flag;flag_case='C';
  elseif length(flag)==1
    if (flag ~= 1), flag_case='A';end
    if (flag == 1), flag_case='B';end
  end
end

if flag_case == 'A'
  be=irf_resamp(b,e);
  
  nl=irf_norm(be); % along the B
  nn=irf_norm(irf_cross(irf_cross(be,v),be)); % closest to given vn vector
  nm=irf_cross(nn,nl); % in (vn x b) direction
  
  % estimate e in new coordinates
  en=irf_dot(e,nn,1);
  el=irf_dot(e,nl,1);
  em=irf_dot(e,nm,1);
  
  %  emp=[e(:,1) el em en];
  emp=e; emp(:,end-2)=el;emp(:,end-1)=em;emp(:,end)=en;
elseif flag_case == 'B'
  nn=irf_norm(v);
  nm=irf_norm(irf_cross(nn,mean(b)));
  nl=irf_cross(nm,nn);
  
  % estimate e in new coordinates
  en=irf_dot(e,nn,1);
  el=irf_dot(e,nl,1);
  em=irf_dot(e,nm,1);
  
  emp=e; emp(:,[end-2 end-1 end])=[el em en];
elseif flag_case == 'C'
  nn=irf_norm(v);
  nm=irf_norm(irf_cross(nn,L_direction));
  nl=irf_cross(nm,nn);
  
  % estimate e in new coordinates
  en=irf_dot(e,nn,1);
  el=irf_dot(e,nl,1);
  em=irf_dot(e,nm,1);
  
  emp=e; emp(:,[end-2 end-1 end])=[el em en];
end
