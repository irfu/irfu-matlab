function [emp,nl,nm,nn]=av_c_mp(e,b,v,flag)
% function       [xmp]=av_c_mp(x,b,v,flag)
% function [xmp,L,M,N]=av_c_mp(x,b,v,flag)
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

if nargin ==3, flag_case='A';end
if (nargin ==4)
  if length(flag)==3,
    L_direction=flag;clear flag;flag_case='C';
  elseif length(flag)==1,
    if (flag ~= 1), flag_case='A';end
    if (flag == 1), flag_case='B';end
  end
end

if flag_case == 'A',
  be=av_interp(b,e);

  nl=av_norm(be); % along the B
  nn=av_norm(av_cross(av_cross(be,v),be)); % closest to given vn vector
  nm=av_cross(nn,nl); % in (vn x b) direction

  % estimate e in new coordinates
  en=av_dot(e,nn,1);
  el=av_dot(e,nl,1);
  em=av_dot(e,nm,1);

%  emp=[e(:,1) el em en];
    emp=e; emp(:,end-2)=el;emp(:,end-1)=em;emp(:,end)=en;
elseif flag_case == 'B',
  nn=av_norm(v);
  nm=av_norm(av_cross(nn,mean(b)));
  nl=av_cross(nm,nn);

  % estimate e in new coordinates
  en=av_dot(e,nn,1);
  el=av_dot(e,nl,1);
  em=av_dot(e,nm,1);

  emp=e; emp(:,[end-2 end-1 end])=[el em en];
elseif flag_case == 'C',
  nn=av_norm(v);
  nm=av_norm(av_cross(nn,L_direction));
  nl=av_cross(nm,nn);

  % estimate e in new coordinates
  en=av_dot(e,nn,1);
  el=av_dot(e,nl,1);
  em=av_dot(e,nm,1);

  emp=e; emp(:,[end-2 end-1 end])=[el em en];
end
