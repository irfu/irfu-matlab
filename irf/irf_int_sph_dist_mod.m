function [pst] = irf_int_sph_dist_mod(F,v,phi,th,vg,varargin)
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
zphat = [0,0,1]; % integrate along this axes in 2D, has no use in 1D
nMC = 10; % number of Monte Carlo iterations
vzint = [-inf,inf]; % limit on out-of-plane velocity
aint = [-180,180]; % limit on out-of-plane velocity
projDim = 1; % number of dimensions of the projection
weight = 'none'; % how number of MC points is weighted to data

args = varargin;
nargs = length(varargin);

% loop to check for flags
have_options = nargs > 1;
while have_options
    switch(lower(args{1}))
        case 'x'
            xphat = args{2};
        case 'z'
            zphat = args{2};
        case 'phig' % azimuthal angle of projection plane
            % If this flag is given, it assumes 2D projection
            phig = args{2};
            projDim = 2;
        case 'nmc'
            nMC = args{2};
        case 'vzint'
            vzint = args{2};
        case 'aint'
            aint = args{2};
        case 'weight'
            weight = args{2};
    end
    args = args(3:end);
    if isempty(args), break, end
end


%% Initiate initiate various things

% complete RH system
yphat = cross(zphat,xphat);

% diffs
dV = diff(v); dV = [dV(1),dV]; % quick and dirty
dPhi = abs(median(diff(phi))); % constant
dTh = abs(median(diff(th))); % constant

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
    nAzg = 1; % for 1D
end

% bin edges
% vg_edges = [vg(1)-dVg(1)/2,vg+dVg/2]; % quick and dirty
vg_edges = zeros(1,length(vg)+1);
vg_edges(1) = vg(1)-diff(vg(1:2))/2;
vg_edges(2:end-1) = vg(1:end-1)+diff(vg)/2;
vg_edges(end) = vg(end)+diff(vg(end-1:end))/2;

if projDim == 2
    phig_edges = [phig-dPhig/2,phig(end)+dPhig/2];
    phig_center = phig;
end

% primed (grid) diffs
dVg = diff(vg_edges);

% convert to cartesian mesh, only for output
if projDim == 2
%     [phiMesh,vMesh] = meshgrid(phig_edges+dPhig/2,vg); % Creates the mesh
%     [vxMesh,vyMesh] = pol2cart(phiMesh-pi/nAzg,vMesh);    % Converts to cartesian
    [phiMesh_edges,vMesh_edges] = meshgrid(phig_edges,vg_edges); % Creates the mesh, edges of bins
    [phiMesh,vMesh] = meshgrid(phig_center,vg); % Creates the mesh, center of bins
    [vxMesh_edges,vyMesh_edges] = pol2cart(phiMesh_edges,vMesh_edges);    % Converts to cartesian, edges
    [vxMesh,vyMesh] = pol2cart(phiMesh-pi/nAzg,vMesh);    % Converts to cartesian, centers
end

% Number of instrument bins
nV = length(v);
nAz = length(phi);
nEle = length(th);


% 3D matrices for instrumental bin centers
TH = repmat(th,nV,1,nAz);       % [phi,th,v]
TH = permute(TH,[1,3,2]);       % [v,phi,th]
PHI = repmat(phi,nV,1,nEle);    % [v,phi,th]
VEL = repmat(v,nAz,1,nEle);     % [phi,v,th]
VEL = permute(VEL,[2,1,3]);     % [v,phi,th]
DV = repmat(dV,nAz,1,nEle);     % [phi,v,th]
DV = permute(DV,[2,1,3]);       % [v,phi,th]

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

% init Fp
%Fg = zeros(nAzg+1,nVg);
Fg = zeros(nAzg,nVg);
% Volume element
dtau = ( VEL.^2.*cos(TH).*DV*dPhi*dTh );
% Area or line element (primed)
dAg = vg.^(projDim-1).*dVg*dPhig;


%% Perform projection
% Loop through all instrument bins
for i = 1:nV % velocity (energy)
    for j = 1:nAz % phi
        for k = 1:nEle % theta
            % generate MC points
            nMCt = NMC(i,j,k); % temporary number
            % Ignore bin if value of F is zero to save computations
            if F(i,j,k) == 0
                continue;
            end
            % first is not random
            % why not random?
            dV_MC = [0;(rand(nMCt-1,1)-.5)*dV(i)];
            dPHI_MC = [0;(rand(nMCt-1,1)-.5)*dPhi];
            dTH_MC = [0;(rand(nMCt-1,1)-.5)*dTh];
            
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
            else
                % convert to polar coordinates (phip could become negative)
                [phip,vp] = cart2pol(vxp,vyp);
                % fix if negative
                phip(phip<0) = 2*pi+phip(phip<0);
            end
            
            % get indicies for all MC points
            iVg = discretize(vp,vg_edges);
            
            if projDim == 2
              % phig - projection grid centers
              % phig_edges - projection grid edges
              % phip - MC particle angle
              % discretize is similar to histcounts, or histcn
              % iAzg can have values from 1 to numel(phi_edges)-1
              iAzg = discretize(phip,phig_edges); 
            else
                iAzg = ones(1,nMCt);
            end
            
            % Loop through MC points and add value of instrument bin to the
            % appropriate projection bin
            for l = 1:nMCt
                % add value to appropriate projection bin
                if usePoint(l) && ~isempty(iAzg(l)) && ~isempty(iVg(l)) && (iAzg(l)<nAzg+1 || iAzg(l)==1) && iVg(l)<nVg+1
                    Fg(iAzg(l),iVg(l)) = Fg(iAzg(l),iVg(l))+F(i,j,k)*dtau(i,j,k)/dAg(iVg(l))/nMCt;
                end
            end
        end
    end
end


%% Output
if projDim == 2
    % fix for interp shading, otherwise the last row can be whatever
    Fg(end,:) = mean([Fg(end-1,:);Fg(1,:)]);
end

% Calculate density
if projDim == 1
    dens = sum(Fg.*dAg);
else
    dAG = repmat(dAg,nAzg,1);
    dens = sum(sum(Fg(1:end,:).*dAG));
end

% Calculate velocity moment

if projDim == 1
    vel = nansum(Fg.*dAg.*vg);
else
    VG = repmat(vg',1,nAzg);
    PHIG = repmat(phig,nVg,1);
    [VXG,VYG] = pol2cart(PHIG,VG);
    
    vel = [0,0];
    for l = 1:nVg
        for m = 1:nAzg
            vel = vel+[VXG(l,m),VYG(l,m)]*Fg(m,l)*dAg(l);
        end
    end
end
vel = vel/dens;


% output
pst = [];
pst.F = Fg;
if projDim == 1
    pst.v = vg;
else
    pst.vx = vxMesh;
    pst.vy = vyMesh;
    pst.vx_edges = vxMesh_edges;
    pst.vy_edges = vyMesh_edges;
end
pst.dens = dens;
pst.vel = vel;

end