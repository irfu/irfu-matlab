function [x, y, z, correlation, intEdt, Bz, B0, dEk, dEn, Ek, En] = irf_match_phibe_dir(B,E,angles,f)
% IRF_MATCH_PHIBE_DIR Get propagation direction by matching dBpar and "phi".
%   Tries different propagation directions and finds the direction
%   perpendicular to the magnetic field that gives the best correlation
%   between the electrostatic potential and the parallel wave magnetic field
%   according to int(E)dt = waveB*B0/n*e*mu0.
%
%   [x,y,z,correlation,intEdt,Bz,dEk,dEn,Ek,En] = IRF_MATCH_PHIBE_DIR(B,E,angles,f)
%
%   Input
%       B - magnetic field (to be filtered if f is given)
%       E - electric field (to be filtered if f is given)
%       angles - the angles in degrees to try (1-180 default)
%       f - filter frequency
%
%   Output
%       x - normal direction (size: n_tries x 3)
%       y - propagation direction
%       z - magnetic field direction.
%       correlation - correlation vector
%       intEdt - "potential"
%       Bz - wave magnetic field in parallel direction
%       B0 - mean magnetic field
%       dEk - wave electric field in propagation direction
%       dEn - wave electric field in propagation normal direction
%       Ek - electric field in propagation direction
%       En - electric field in propagation normal direction
%
%   Examples:
%       % Direction
%       angles=1:3:360;
%       f_highpass=7;
%       [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=IRF_MATCH_PHIBE_DIR(B,E,angles,f_highpass);
%       i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
%       direction=x(i_dir,:);
%
%       % Velocity and density
%       n=linspace(0.01,0.1,100);
%       v=linspace(200,2000,100);
%       [corr_v,phi_E,phi_B]=irf_match_phibe_v(B0,Bz,intEdt(:,[1 1+i_dir]),n,v);
%       i_v=find(corr_v(:,1)==min(corr_v(:,1)));
%       velocity=v(i_v);
%
%       % Figures
%       gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En);
%       imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);
%
%       i_n=50; % if more than one densitiy, choose one by specifying index
%       gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
%       imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);
%
%       figure; h=axes;
%       axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);
%
%   See also IRF_MATCH_PHIBE_V, IRF_MATCH_PHIBE_VIS

% Resample B to E if they have different size
if size(B,1) ~= size(E,1)
  B=irf_resamp(B,E);
end

% Filter if f is given, otherwise assume it is filtered
if exist('f','var')
  BAC=irf_filt(B,f,0,450,5);
  EAC=irf_filt(E,f,0,450,5);
else
  BAC=B;
  EAC=E;
end

% Get background magnetic field, for irf_match_phibe_v
if size(B,2)==4; B0=mean(irf_abs(B,1));
elseif size(B,2)==5; B0=mean(B(:,end));
end

% If no angles are specified, set 1,4,7,...,158 as default
if ~exist('angles','var')
  angles=1:3:360;
end
na=length(angles); % number of angles

% Set up coordinate systems
if 1
  z=repmat(irf_norm(mean(B(:,2:4),1)),na,1); % B/z direction, tries*3
  y=irf_norm(irf_cross(irf_cross(z,[1 0 0]),z)); % perp1
  x=irf_norm(irf_cross(y,z)); % perp2

  theta=(0:2*pi/na:2*pi-pi/na)'; % angles
  xn=irf_norm(x.*repmat(cos(theta),1,3)+y.*repmat(sin(theta),1,3));
  y=irf_cross(z,xn);
  x=xn;
end

% Field aligned B
Bz=irf_dot(BAC,z(1,:));

% Allocate correlations
correlation=zeros(na,1);

% Allocate vectors, 4 first used for illustration
dEk=[EAC(:,1) zeros(size(E,1),na)]; % field to integrate
dEn=[EAC(:,1) zeros(size(E,1),na)]; % normal field
Ek=[E(:,1) zeros(size(E,1),na)]; % unfiltered
En=[E(:,1) zeros(size(E,1),na)]; % unfiltered
intEdt=[EAC(:,1) zeros(size(E,1),na)]; % potential

% Integrate E in all x-directions
for k=1:na
  % Used for visualization
  dEk(:,k+1)=irf_dot(EAC(:,(2:4)),x(k,:));
  dEn(:,k+1)=irf_dot(EAC(:,(2:4)),y(k,:));
  Ek(:,k+1)=irf_dot(E(:,(2:4)),x(k,:));
  En(:,k+1)=irf_dot(E(:,(2:4)),y(k,:));

  % Get Phi_E = int(Ek), there's no minus since the field is integrated
  % in the opposite direction of the wave propagation direction.
  prel=irf_integrate([dEk(:,1) dEk(:,k+1)]);
  intEdt(:,k+1)=prel(:,2)-mean(prel(:,2));

  % Get correlation
  correlation(k,1)=xcorr(intEdt(:,k+1),Bz(:,2),0,'coeff');
end


if 0
  vis.type='direction'; %#ok<UNRCH>
  vis.x=x;
  vis.y=y;
  vis.z=z;
  vis.correlation=correlation;
  vis.intEdt=intEdt;
  vis.Ek=Ek;
  vis.En=En;
  vis.ufEk=ufEk;
  vis.ufEn=ufEn;
end