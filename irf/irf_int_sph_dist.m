function [pst] = irf_int_sph_dist(F,v,phi,th,vg,varargin)
% IRF_INT_SPH_DIST Integrate a spherical distribution function to a line/plane.
%
%   pst = IRF_INT_SPH_DIST(F,v,phi,vg,'Opt1',OptVal1,...) Calculates
%   the integrated distribution function of a 3D distribution function in a
%   spherical grid. F is a 3D matrix containing values of the flux (index
%   order of F is [velocity,phi,th]), phi is the azimuthal and th the
%   elevation angles in radians, and vg is the bin centers of the velocity
%   of the projection grid.
%
%   Options:
%   'phig'  -   the bin centers of the azimuthal angle of the projection in
%               radians in the span [0,2*pi]. If this input is given, the
%               projection will be 2D. If it is omitted, the projection
%               will be 1D.
%   'x'     -   axis that is not integrated along in 1D and x-axis
%               (phig = 0) in 2D. x = [1,0,0] if omitted.
%   'z'     -   Axis that is integrated along in 2D. z = [0,0,1] if omitted.
%   'nMC'   -   average number of Monte Carlo iterations used per bin for
%               integration, default is 10. Number of iterations can be
%               weighted data value in each bin.
%   'weight'-   how the number of MC iterations per bin is weighted, can be
%               'none' (default), 'lin' or 'log'
%   'vzint' -   set limits on the out-of-plane velocity interval in 2D and
%               "transverse" velocity in 1D.
%   'aint'  -   angular limit in degrees, can be combined with vzlim
%   'base'  -   coordinate base, 'pol' for polar or 'cart' for cartesian,
%               does not matter in 1D case
%   've'    -   velocity edges of instrument (from delta_energy) with same
%               units as v and one element longer than v
%   'dphi'  -   array of dphi in cases where phi are not evenly separated
%   'dth'   -   array of dth in cases where phi are not evenly separated
%   'counts'-   if F is a PDist converted to counts, this flag should be
%               set to 1. It is 0 by default.
%
%
%   Output is a structure that contains the fields:
%   'F'     -   integrated flux. For 2D, F(end,:) contains no information
%               but is there so it can be plotted.
%   'v'     -   if 1D, velocity of grid centers (same as vg).
%   'vx'/vy -   if 2D, velocity of grid centers.
%   'dens'  -   zeroth order moment of projected distribution (number
%               density).
%   'vel'   -   first order moment of projected distribution divided by the
%               zeroth order moment (flow velocity).
%
%   Description of method:
%   The function goes through all instrument bins and finds the best match
%   on the projection bin. The value in the projection bin, F, is updated
%   as F = F+f*dTau/dA, where f is the instrument value in the bin, dTau is
%   the volume element in velocity space of the instrument bin and dA is
%   the area or line element in the projection bin.
%   An instrument bin can actually cover several projection bin and the
%   value should be added to all those bins, scaled to the area the
%   instrument bin covers of a given projection bin. This area is
%   calculated with a Monte Carlo integration scheme, where many
%   "particles" are generated somewhere inside the instrument bin, each
%   assigned a fraction of f. The "particles" are then projected onto the
%   line or plane.
%
%   See also: MMS.PLOT_INT_PROJECTION

%   Written by: Andreas Johlander, andreasj@irfu.se
%
%   TODO:   Optimization


%% Check for flags in input
% Set default values
xphat = [1,0,0]; % axes projection is done against in 1D, x-axis in 2D
yphat = [0,1,0];
zphat = [0,0,1]; % integrate along this axes in 2D, has no use in 1D
nMC = 10; % number of Monte Carlo iterations
vzint = [-inf,inf]; % limit on out-of-plane velocity
aint = [-180,180]; % limit on out-of-plane velocity
projDim = 1; % number of dimensions of the projection
weight = 'none'; % how number of MC points is weighted to data
base = 'pol'; % If 1D then this does not matter
veInput = 0; % input energy differences
veInputEdges = 0; %
counts = 0;

args = varargin;
nargs = length(varargin);

dPhi = abs(median(diff(phi)))*ones(size(phi)); % constant
dTh = abs(median(diff(th)))*ones(size(th)); % constant

% loop to check for flags
have_options = nargs > 1;
while have_options
  switch(lower(args{1}))
    case 'x'
      xphat = args{2};
    case 'z'
      zphat = args{2};
      % If this flag is given, it assumes 2D projection
      projDim = 2;
    case 'phig' % azimuthal angle of projection plane
      phig = args{2};
    case 'nmc'
      nMC = args{2};
    case 'vzint'
      vzint = args{2};
    case 'aint'
      aint = args{2};
    case 'weight'
      weight = args{2};
    case 'base'
      base = args{2};
    case 've'
      ve = args{2};
      veInput = 1;
    case 'vg_edges'
      vg_edges = args{2};
      veInputEdges = 1;
    case 'dphi'
      dPhi = args{2};
    case 'dth'
      dTh = args{2};
    case 'counts'
      counts = args{2};
  end
  args = args(3:end);
  if isempty(args), break, end
end

%% Initiate initiate various things

% complete RH system
xphat = xphat./sqrt(sum(xphat.^2));
if ~isequal(xphat,zphat)
  yphat = cross(zphat,xphat); % zphat define as default [0 0 1] or read in as optional input above
  yphat = yphat./sqrt(sum(yphat.^2));
  zphat = cross(xphat,yphat); % z = cross(x,cross(z,x)) % enforce z to be orthogonal to x
  zphat = zphat./sqrt(sum(zphat.^2));
else
  zphat = cross(xphat,yphat);
  zphat = zphat./sqrt(sum(zphat.^2));
  yphat = cross(zphat,xphat);
  yphat = yphat./sqrt(sum(yphat.^2));
end

% diffs of instrument bins
% velocity
if veInput
  dVm = v-ve(1:end-1); dVp = ve(2:end)-v; % % minus and plus velocity from center
  dV = dVm+dVp; % total difference
else
  dV = diff(v); dV = [dV(1),dV]; % quick and dirty
  dVm = diff(v)/2; % dVp = diff(v)/2; % minus and plus velocity from center
end

% primed (grid) diffs
% dVg = diff(vg); dVg = [dVg(1),dVg]; % quick and dirty
if projDim == 2
  dPhig = median(diff(phig)); % constant
else
  dPhig = 1; % unity for 1D
end
% Number of projection bins
nVg = length(vg);
if projDim == 2
  nAzg = length(phig);
else
  nAzg = 0; % for 1D
end

% Number of instrument bins
nAz = length(phi);
nEle = length(th);
nV = length(v);



%% bin edges
if veInputEdges % get vg from vg_edges
  vg_edges = vg_edges;
  vg = vg_edges(1:(end-1)) + 0.5*diff(vg_edges);
  nVg = numel(vg);
else % get vg_edges from vg
  vg_edges = zeros(1,length(vg)+1);
  vg_edges(1) = vg(1)-diff(vg(1:2))/2;
  vg_edges(2:end-1) = vg(1:end-1)+diff(vg)/2;
  vg_edges(end) = vg(end)+diff(vg(end-1:end))/2;
end

switch lower(base)
  case 'pol'
    if projDim == 2
      phig_edges = [phig-dPhig/2,phig(end)+dPhig/2];
    end

    % primed (grid) diffs
    dVg = diff(vg_edges);

    % convert to cartesian mesh, only for output
    if projDim == 2
      [phiMesh,vMesh] = meshgrid(phig_edges+dPhig/2,vg); % Creates the mesh, center of bins, phi has one extra bin at the end
      [vxMesh,vyMesh] = pol2cart(phiMesh-pi/nAzg,vMesh);    % Converts to cartesian

      [phiMesh_edges,vMesh_edges] = meshgrid(phig_edges,vg_edges); % Creates the mesh, edges of bins
      [vxMesh_edges,vyMesh_edges] = pol2cart(phiMesh_edges,vMesh_edges); % Converts to cartesian, edges
    end



  case 'cart'
    % for cartesian grid, the velocity bins must all be equal
    % a linearly spaced grid can have small roundoff differences in step
    % with std, there could potentially be some outlier? i dont know
    meandiff = mean(diff(vg));
    errtol = 1e-2; % 1%
    if not(all((diff(vg)/meandiff-1)<errtol))
      error('For a cartesian grid (default), all velocity bins diff(vg) must be equal.');
    end
    dVg = vg(2)-vg(1);
end

% 3D matrices for instrumental bin centers
TH = repmat(th,nV,1,nAz);       % [phi,th,v]
TH = permute(TH,[1,3,2]);       % [v,phi,th]
PHI = repmat(phi,nV,1,nEle);    % [v,phi,th]
VEL = repmat(v,nAz,1,nEle);     % [phi,v,th]
VEL = permute(VEL,[2,1,3]);     % [v,phi,th]
DV = repmat(dV,nAz,1,nEle);     % [phi,v,th]
DV = permute(DV,[2,1,3]);       % [v,phi,th]
DTH = repmat(dTh,nV,1,nAz);       % [phi,th,v]
DTH = permute(DTH,[1,3,2]);       % [v,phi,th]
DPHI = repmat(dPhi,nV,1,nEle);   % [v,phi,th]

% Weighting of number of Monte Carlo particles
Nsum = nMC*numel(find(F)); % total number of Monte Carlo particles
switch weight
  case 'none'
    % 3D matrix with values of nMC for each bin
    NMC = zeros(size(F)); % no points when data is 0
    NMC(F~=0) = nMC;
  case 'lin'
    NMC = ceil(Nsum/sum(sum(sum(F)))*F);
  case 'log'
    NMC = ceil(Nsum/(sum(sum(sum(log10(F+1)))))*log10(F+1));
end



% "Volume" element in velocity space of an instrument bin
dtau = ( VEL.^2.*cos(TH).*DV.*DPHI.*DTH );
% set grid data matrix and grid "area" element
switch lower(base)
  case 'pol'
    % init Fp
    Fg = zeros(nAzg+1,nVg);
    Fg_ = zeros(nAzg,nVg); % use this one with 'edges bins'
    % Area or line element (primed)
    dAg = vg.^(projDim-1).*dVg*dPhig;

  case 'cart'
    Fg = zeros(nVg,nVg);
    dAg = dVg^2;
end

%% Perform projection
% Loop through all instrument bins
for i = 1:nV % velocity (energy)
  for j = 1:nAz % phi
    for k = 1:nEle % theta
      % generate MC points
      nMCt = NMC(i,j,k); % temporary number
      % Ignore bin if value of F is zero to save computations
      if F(i,j,k) == 0 || isnan(F(i,j,k))
        continue;
      end

      % Construct Monte Carlo particles within particle bins
      % first is not random
      dV_MC = [0;-rand(nMCt-1,1)*dV(i)-dVm(1)]; % velocity within [-dVm,+dVp]
      dPHI_MC = [0;(rand(nMCt-1,1)-.5)*dPhi(j)];
      dTH_MC = [0;(rand(nMCt-1,1)-.5)*dTh(k)];

      % convert instrument bin to cartesian velocity
      [vx,vy,vz] = sph2cart(PHI(i,j,k)+dPHI_MC,TH(i,j,k)+dTH_MC,VEL(i,j,k)+dV_MC);

      % Get velocities in primed coordinate system
      vxp = [vx,vy,vz]*xphat'; % all MC points
      vyp = [vx,vy,vz]*yphat';
      vzp = [vx,vy,vz]*zphat'; % all MC points
      vabsp = sqrt(vxp.^2+vyp.^2+vzp.^2);
      if projDim == 1 % get transverse velocity sqrt(vy^2+vz^2)
        vzp = sqrt(vyp.^2+vzp.^2); % call it vzp
      end
      alpha = asind(vzp./vabsp);

      % If "particle" is outside allowed interval, don't use point
      usePoint = (vzp >= vzint(1) & vzp <= vzint(2) & alpha >= aint(1) & alpha <= aint(2));

      if projDim == 1
        vp = vxp;
      elseif strcmp(base,'pol')
        % convert to polar coordinates (phip could become negative)
        [phip,vp] = cart2pol(vxp,vyp);
        % fix if negative
        phip(phip<0) = 2*pi+phip(phip<0);
      end


      % different procedure for 1D or polar OR cartesian
      if strcmpi(base,'pol') || projDim == 1
        % ------ 1D AND POLAR CASE ------
        % get indicies for all MC points
        iVg = discretize(vp,vg_edges);
        % fixes bug that exists on some systems, may influence
        % performance
        iVg(iVg==0) = nan;

        if projDim == 2
          iAzg = discretize(phip,phig_edges);
        else
          iAzg = ones(1,nMCt);
        end

        % Loop through MC points and add value of instrument bin to the
        % appropriate projection bin
        for l = 1:nMCt
          % add value to appropriate projection bin
          if usePoint(l) && ~isempty(iAzg(l)) && ~isempty(iVg(l)) && (iAzg(l)<nAzg+1 || iAzg(l)==1) && iVg(l)<nVg+1
            if counts
              Fg(iAzg(l),iVg(l)) = Fg(iAzg(l),iVg(l))+F(i,j,k)/nMCt; % no dtau/dAg weighting
            else
              Fg(iAzg(l),iVg(l)) = Fg(iAzg(l),iVg(l))+F(i,j,k)*dtau(i,j,k)/dAg(iVg(l))/nMCt;
            end
            if projDim == 2
              if counts
                Fg_(iAzg(l),iVg(l)) = Fg(iAzg(l),iVg(l))+F(i,j,k)/nMCt; % no dtau/dAg weighting
              else
                Fg_(iAzg(l),iVg(l)) = Fg(iAzg(l),iVg(l))+F(i,j,k)*dtau(i,j,k)/dAg(iVg(l))/nMCt;
              end
            end
          end
        end

      elseif strcmpi(base,'cart')
        % ------ CARTESIAN CASE ------

        % get indicies for all MC points
        iVxg = discretize(vxp,vg_edges);
        iVyg = discretize(vyp,vg_edges);
        % fixes bug that exists on some systems, may influence
        % performance
        iVxg(iVxg==0) = nan;
        iVyg(iVyg==0) = nan;

        % Loop through MC points and add value of instrument bin to the
        % appropriate projection bin
        for l = 1:nMCt
          if usePoint(l) && vxp(l)>min(vg_edges) && vxp(l)<max(vg_edges) && vyp(l)>min(vg_edges) && vyp(l)<max(vg_edges)
            if counts
              Fg(iVxg(l),iVyg(l)) = Fg(iVxg(l),iVyg(l))+F(i,j,k)/nMCt; % no dtau/dAg weighting
            else
              Fg(iVxg(l),iVyg(l)) = Fg(iVxg(l),iVyg(l))+F(i,j,k)*dtau(i,j,k)/dAg/nMCt;
            end
          end
        end
      end

    end
  end
end


%% Output
if projDim == 2 && strcmpi(base,'pol')
  % fix for interp shading, otherwise the last row can be whatever
  Fg(end,:) = mean([Fg(end-1,:);Fg(1,:)]);
end

% Calculate density
if projDim == 1
  dens = sum(Fg.*dAg);
elseif strcmpi(base,'pol')
  dAG = repmat(dAg,nAzg,1);
  dens = sum(sum(Fg(1:end-1,:).*dAG));
elseif strcmpi(base,'cart')
  dens = sum(sum(Fg*dAg));
end

% Calculate velocity moment

if projDim == 1
  vel = sum(Fg.*dAg.*vg,'omitnan');
elseif strcmpi(base,'pol')
  VG = repmat(vg',1,nAzg);
  PHIG = repmat(phig,nVg,1);
  [VXG,VYG] = pol2cart(PHIG,VG);

  vel = [0,0];
  for l = 1:nVg
    for m = 1:nAzg
      vel = vel+[VXG(l,m),VYG(l,m)]*Fg(m,l)*dAg(l);
    end
  end
elseif strcmpi(base,'cart')
  vel = [0,0]; % whatever
end
vel = vel/dens;


% output
pst = [];
pst.F = Fg;
if projDim == 1
  pst.v = vg;
  pst.v_edges = vg_edges;
elseif strcmpi(base,'pol')
  pst.vx = vxMesh;
  pst.vy = vyMesh;
  pst.F_using_edges = Fg_;
  pst.vx_edges = vxMesh_edges;
  pst.vy_edges = vyMesh_edges;
elseif strcmpi(base,'cart')
  pst.vx = vg;
  pst.vy = vg;
  pst.F_using_edges = [[Fg,zeros(nVg,1)];zeros(1,nVg+1)];
  pst.vx_edges = vg_edges;
  pst.vy_edges = vg_edges;

end
pst.dens = dens;
pst.vel = vel;

end