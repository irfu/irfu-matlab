function [jz,jp,nj,angle]=irf_jz(v,B,dB,deg_p,deg_z,n_av)
%IRF_JZ   Estimate the current given velocity and magnetic field
%
% [jz,jp,nj,angle]=irf_jz(v,B,dB,deg_p,deg_z,n_av);        % method A
% [jz,jp,nj,angle]=irf_jz(v,B);                            % method A
% [j]=irf_jz(v,B,dB,n_av);                                 % method B
% [j]=irf_jz(v,B);                                         % method B
% [jz,jp]=irf_jz(v,B);                                     % method B
%
% Estimates the current given velocity and magnetic field.
% Two methods are implemented:
% A) assume that spacecraft are crossing plane current sheets (tangential
%    discontinuities) that can have normals at an angle with respect to the
%    velocity. The normal is found from cross-product of magnetic field
%    changes and the background magnetic field direction.
%    This method is useful for example for auroral current sheet crossings.
% B) neglect changes in B along the velocity vector
%    (thus all current is perpendicular to velocity). Useful for example for MP crossings.
%     j=(1/mu0/v^2/dt)   v x dtB
%
% v     [[t] vx vy vz] velocity of the spacecraft [km/s] or
%       when medium velocity dominates the negative of medium velocity
% B     [t Bx By Bz] - background magnetic field [nT]
% dB    [t dBx dBy dBz] - disturbance magnetic field [nT] (if not given uses background)
% n_av - how many points of B to average for current estimates (in practice
%        this is lowpass filter with filter frequency at Fs/2/n_av)
%
%        for case 'A' angle between current sheet and velocity should be at least deg_p degrees
%       (in a plane perp to B) and v pitch angle should be at least deg_z,
%        otherwise current is set to NaN  [default 10 degrees]
%
% j  - current [A/m2] total current vector.
% jz - current [A/m2] parallel to B, time from B, pozitive when jz along B
% jp - current [A/m2] perpendicular to B, time from B, pozitive when |B| increases along spacecraft trajectory
% nj - the vector direction of the current sheet
% angle - azimuthal angle between current sheet and velocity component perpendicular to B
%         angle is measured anticlockwise in XY plane from X axis
%         system is defined such that B is Z and V_perp is X
%
% See also C_PL_J, C_4_J

if size(v,2)==3, v=[v(:,1)*0+B(1,1) v];end

if nargout<=2, method='B'; else, method='A';end
if nargin <3, dB=B;end
flag_average=0;
if strcmp(method,'A')
  if nargin==6, flag_average=1;end
  if nargin <5, deg_z=10;end
  if nargin <4, deg_p=10;end
elseif strcmp(method,'B')
  if nargin==4, flag_average=1;n_av=deg_p;end
end


if flag_average == 1
  dBf=irf_filt(dB,0,1/n_av,2,3);
  tt=dBf(:,1);
  [~ ,dtB] = gradient(dBf(:,2:end),1,tt); 
  dtB      = [tt dtB];
else
  tt=dB(:,1);
  % create dtB
  [~ ,dtB] = gradient(dB(:,2:end),1,tt); 
  dtB      = [tt dtB];
end

muo=4*pi/1e7;

switch method
  case 'A'
    nb=irf_resamp(irf_norm(B),dtB);
    vb=irf_resamp(v,dtB);
    dtBnj=irf_cross(dtB,nb);
    dtBnb=irf_dot(dtB,nb);
    nj=irf_norm(dtBnj); % current sheet plane normal

    vnj=irf_dot(vb,nj); % v component in the direction normal to current sheet
    vp=irf_cross(nb,irf_cross(vb,nb)); % velocity perp to B
    sin_angle=irf_dot(irf_cross(irf_norm(vp),nj),nb);xxx=irf_dot(vp,nj);ind=find(xxx(:,2)<0);sin_angle(ind,2)=-sin_angle(ind,2);
    angle=[dtB(:,1) -asin(sin_angle(:,2))*180/pi];

    absvp=irf_abs(vp,1); absvb=irf_abs(vb,1);absvnj=abs(vnj(:,2));
    indnan_z=find(absvp < absvb*sin(pi/180*deg_z));
    indnan_p=find(absvnj < absvp*sin(pi/180*deg_p));
    jz=[dtB(:,1) irf_abs(dtBnj,1)./vnj(:,2)/muo*1e-12];jz(indnan_z,2)=NaN;jz(indnan_p,2)=NaN;
    jp=[dtB(:,1) abs(dtBnb(:,2))./vnj(:,2)/muo*1e-12];jp(indnan_z,2)=NaN;jp(indnan_p,2)=NaN;
  case 'B'
    dt=dB(2,1)-dB(1,1);
    irf.log('warning',['fs=' num2str(1/dt,3) ' irf_jz()']);
    vb=irf_resamp(v,dtB);
    jxx=irf_vec_x_scal(irf_cross(dtB,vb),[vb(:,1) irf_abs(vb,1)],-2);
    j=irf_tappl(jxx,'*(-1)/(4*pi/1e7)*1e-12');
    jz=irf_dot(j,irf_norm(B));
    jp=irf_dot(j,irf_norm(irf_cross(B,vb)));
    if nargout==1, jz=j; end
end
