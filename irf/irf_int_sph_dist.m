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
%   'z'     -   Axis that is integrated along in 2D. Has no use in 1D.
%               z = [0,0,1] if omitted.
%   'nMC'   -   number of Monte Carlo iterations used for integration,
%               default is 10.
%   'vzint' -   set limits on the out-of-plane velocity interval in 2D and
%               "transverse" velocity in 1D.
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
projDim = 1; % number of dimensions of the projection

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
    end
    args = args(3:end);
    if isempty(args), break, end
end


% complete RH system
yphat = cross(zphat,xphat);

% diffs
dV = diff(v); dV = [dV(1),dV]; % quick and dirty
dPhi = abs(median(diff(phi))); % constant
dTh = abs(median(diff(th))); % constant

% primed (grid) diffs
dVg = diff(vg); dVg = [dVg(1),dVg]; % quick and dirty
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

% bin edges
vg_edges = [vg(1)-dVg(1)/2,vg+dVg/2]; % quick and dirty
if projDim == 2
    phig_edges = [phig-dPhig/2,phig(end)+dPhig/2];
end

% convert to cartesian mesh, only for output
if projDim == 2
    [phiMesh,vMesh] = meshgrid(phig_edges+dPhig/2,vg); % Creates the mesh
    [vxMesh,vyMesh] = pol2cart(phiMesh-pi/nAzg,vMesh);    % Converts to cartesian
end

% Number of instrument bins
nV = length(v);
nAz = length(phi);
nEle = length(th);


% 3D matrices for bin centers
TH = repmat(th,nV,1,nAz);       % [phi,th,v]
TH = permute(TH,[1,3,2]);       % [v,phi,th]
PHI = repmat(phi,nV,1,nEle);    % [v,phi,th]
VEL = repmat(v,nAz,1,nEle);     % [phi,v,th]
VEL = permute(VEL,[2,1,3]);     % [v,phi,th]
DV = repmat(dV,nAz,1,nEle);     % [phi,v,th]
DV = permute(DV,[2,1,3]);       % [v,phi,th]


% init Fp
Fg = zeros(nAzg+1,nVg);
% Volume element
dtau = ( VEL.^2.*cos(TH).*DV*dPhi*dTh );
% Area or line element (primed)
dAg = vg.^(projDim-1).*dVg*dPhig;

% Loop through all instrument bins
for i = 1:nV % velocity (energy)
    for j = 1:nAz % phi
        for k = 1:nEle % theta
            % generate MC points
            % first is not random
            dV_MC = [0;(rand(nMC-1,1)-.5)*dV(i)];
            dPHI_MC = [0;(rand(nMC-1,1)-.5)*dPhi];
            dTH_MC = [0;(rand(nMC-1,1)-.5)*dTh];
            
            % convert instrument bin to cartesian velocity
            [vx,vy,vz] = sph2cart(PHI(i,j,k)+dPHI_MC,TH(i,j,k)+dTH_MC,VEL(i,j,k)+dV_MC);
            
            % Get velocities in primed coordinate system
            vxp = sum([vx,vy,vz].*xphat,2); % all MC points
            if projDim == 1 % get transverse velocity sqrt(vy^2+vz^2)
                vzp = dot([vx(1),vy(1),vz(1)],zphat); % only bin center
                vyp = dot([vx(1),vy(1),vz(1)],yphat); % only bin center
                vzp = sqrt(vyp^2+vzp^2); % call it vzp
            else % get y and z for 2D
                vyp = sum([vx,vy,vz].*yphat,2); % all MC points
                vzp = dot([vx(1),vy(1),vz(1)],zphat); % only bin center
            end
            
            % If bin center outside allowed interval, set F to zero
            if vzp < vzint(1) || vzp > vzint(2); F(i,j,k) = 0; end
            
            if projDim == 1
                vp = vxp;
            else
                % convert to polar coordinates (phip could become negative)
                [phip,vp] = cart2pol(vxp,vyp);
                % fix if negative
                phip(phip<0) = 2*pi+phip(phip<0);
            end
            
            % not so good but better than throwing away data? Should be
            % very rare anyway.
            
            % Loop through MC points and add value of instrument bin to the
            % appropriate projection bin
            for l = 1:nMC
                iVg = find(vp(l)>vg_edges,1,'last');
                % Add to closest bin if it falls outside
                if isempty(iVg) && vp(l)<vg_edges(1); iVg = 1; end
                if iVg == nVg+1 && vp(l)>vg_edges(end); iVg = nVg; end
                    
                if projDim == 2
                    iAzg = find(phip(l)>phig_edges,1,'last');
                else
                    iAzg = 1;
                end
                
                try % add value to appropriate projection bin
                    Fg(iAzg,iVg) = Fg(iAzg,iVg)+F(i,j,k)*dtau(i,j,k)/dAg(iVg)/nMC;
                catch
                    irf.log('w',['Something went wrong. iAzg = ',num2str(iAzg),',  iVg = ',num2str(iVg)])
                end
            end
        end
    end
end

if projDim == 2
    % fix for interp shading, otherwise the last row can be whatever
    Fg(end,:) = mean([Fg(end-1,:);Fg(1,:)]);
end

% Calculate density
if projDim == 1
    dens = nansum(Fg.*dAg);
else
    dAG = repmat(dAg,nAzg,1);
    dens = nansum(nansum(Fg(1:end-1,:).*dAG));
end

% Calculate velocity moment

if projDim == 1
    vel = nansum(Fg.*dAg.*vg);
else
    VG = repmat(vg',1,nVg);
    PHIG = repmat(phig,nAzg,1);
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
end
pst.dens = dens;
pst.vel = vel;

end