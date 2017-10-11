function [ax,dens,velmom,plspec] = plot_int_projection(varargin)
%MMS.PLOT_INT_PROJECTION Plots integrated FPI data on a specified plane.
%
%   MMS.PLOT_INT_PROJECTION(dist,'Opt1',OptVal1,...) plots projection of
%   particle data in the PDist object dist.
%
%   MMS.PLOT_INT_PROJECTION(AX,...) plots in axes AX instead of current
%   axes.
%
%   [ax,dens,velmom,plspec] = MMS.PLOT_INT_PROJECTION(...) returns axis
%   handle, density, 2D velocity derived from the projected data in SI
%   units and a structure plspec with data used for plotting.
%   plspec contains:
%       F   - the particle flux, for PSD, typically the unit is [s^2m^-5],
%           which is the same as [s^2cm^-3^km^-2]
%       vx  - vx grid in [m/s]
%       vy  - vy grid in [m/s]
%   To recreate plot elsewhere use: 
%       pcolor(plspec.vx,plspec.vy,log10(plspec.F'))
% 
%   Options:
%   'tint'/'t'  - plots data for closest time given an EpochTT object. Time
%               interval not implemented
%   'xyz'       - 3x3 matrix with [x;y;z]. z is normal to the plotted plane and
%               x and y are made orthogonal to z and each other if they are 
%               not already. If you want to plot different planes you have to
%               rotate this matrix -> [y;z;x] -> [z;x;y]. Make sure the system
%               is right-handed.
%   'clim'      - [cmin cmax], colorbar limits in logscale
%   'vlim'      - vmax in [km/s], zoom in to XLim = YLim = vlim*[-1 1]
%   'vzint'     - set limits on the out-of-plane velocity interval [km/s],
%               useful for i.e. electron crescent distributions
%   'nMC'       - number of Monte Carlo iterations used for integration,
%               default is 100
%   'vlabel'    - 1x3 cell array containing strings for axis labels
%               corresponding to x, y, and z
%   'flipx'/'flipy' - boolean value where 1 flips the x/y axis 
%
%   Description of method:
%   The function goes through all instrument bins and finds the best match
%   on the projection bin. The value in the projection bin, F, is updated
%   as F = F+f*dTau/dA, where f is the instrument value in the bin, dTau is
%   the volume element in velocity space of the instrument bin and dA is
%   the area element in the projection bin. 
%   An instrument bin can actually cover several projection bin and the
%   value should be added to all those bins, scaled to the area the
%   instrument bin covers of a given projection bin. This area is
%   calculated with a Monte Carlo integration scheme, where many
%   "particles" are generated somewhere inside the instrument bin, each
%   assigned a fraction of f. The "particles" are then projected onto the
%   plane.

%   Written by: Andreas Johlander, andreasj@irfu.se
% 
%   TODO:   Add time interval averaging
%           Recalculate energies given spacecraft potential input
%           In great need of optimization
%
%   Generally, the code deals with quantities in SI units until plotting.

%% Input

[ax,args,nargs] = axescheck(varargin{:});

irf.log('warning','Please verify that you think the projection is done properly!');

dist = args{1};
args = args(2:end);

% Check if it's electron or ions
if strfind(dist.name,'des')
  isDes = 1; 
elseif strfind(dist.name,'dis')
  isDes = 0;
else
  irf.log('warning','Can''t recognize if input is electron or ions. Assuming it''s electrons.');
  isDes = 1;
end

if isempty(dist); irf.log('warning','Empty input.'); return; end

%Check if data is skymap
if ~strcmp(dist.type,'skymap')
    irf.log('critical','PDist must be skymap format.');
    return;
end


%% Check for flags

% Set default values
doFlipX = 0;
doFlipY = 0;
have_vlabels = 0;
have_clim = 0;
have_vlim = 0;
nMC = 100; % number of Monte Carlo iterations
vzint = [-inf,inf];

have_options = nargs > 1;
while have_options
    switch(lower(args{1}))
        case 't' % time
            t = args{2};
        case 'xyz' % 3x3 matrix where xphat = coord_sys(1,:),...
            coord_sys = args{2};
            xphat = coord_sys(1,:)/norm(coord_sys(1,:)); % new x-axis (x-primed)
            zphat = coord_sys(3,:)/norm(coord_sys(3,:)); % vector to integrate along (z-primed)
            % only x and z are relvant
            yphat = cross(zphat,xphat); yphat = yphat/norm(yphat);
            if abs(acosd(yphat*(coord_sys(2,:)/norm(coord_sys(2,:)))'))>1
                irf.log('warning',['y (perp1) changed from [' num2str(coord_sys(2,:)/norm(coord_sys(2,:)),'% .2f') '] to [' num2str(yphat,'% .2f') '].']);
            end
        case 'clim' 
            clim = args{2};
            have_clim = 1;
        case 'vlim'
            vlim = args{2};
            have_vlim = 1;
        case 'vlabel'
            have_vlabels = 1;
            vlabels = args{2};
        case 'flipx'
            doFlipX = args{2};
        case 'flipy'
            doFlipY = args{2};
        case 'nmc' % number of Monte Carlo iterations
            nMC = args{2};
        case 'vzint' % limit on out-of-plane velocity
            vzint = args{2};
            
    end
    args = args(3:end);
    if isempty(args), break, end
end

if have_vlabels
    vlabelx = vlabels{1};
    vlabely = vlabels{2};
    vlabelz = vlabels{3};
else
    vlabelx = ['v_{x=[' num2str(xphat,'% .2f') ']}'];
    vlabely = ['v_{y=[' num2str(yphat,'% .2f') ']}'];
    vlabelz = ['v_{z=[' num2str(zphat,'% .2f') ']}'];
end



%% Get angles and velocities for spherical instrument grid

u = irf_units;

if isDes == 1; M = u.me; else; M = u.mp; end

% find closest time
it = anjo.fci(t.epochUnix,dist.time.epochUnix);

% 3d data matrix for time index it
F3d = squeeze(double(dist.data(it,:,:,:)))*1e12; % s^3/m^6


emat = dist.energy; % in eV
energy = emat(it,:);
v = sqrt(2*energy*u.e/M); % m/s

if length(v) ~= 32
    error('something went wrong')
end

% azimuthal angle
phi = dist.depend{2}(it,:); % in degrees
%phi = phi+180;
%phi(phi>360) = phi(phi>360)-360;
phi = phi-180;
phi = phi*pi/180; % in radians

if length(phi) ~= 32
    error('something went wrong')
end

% elevation angle
th = dist.depend{3}; % polar angle in degrees
th = th-90; % elevation angle in degrees
th = th*pi/180; % in radians

if length(th) ~= 16
    error('something went wrong')
end


%% Set projection grid
nAzg = 32;

% diffs
dPhig = 2*pi/nAzg;

% bin centers
phig = linspace(0,2*pi-dPhig,nAzg)+dPhig/2;
vg = v; % same as instrument

%% perform projection
[Fg,vx_mesh,vy_mesh,dens,velmom] = project_sph_distr(F3d,v,phi,th,zphat,xphat,vg,phig,nMC,vzint*1e3);

% put nans instead of 0s
Fg(Fg==0) = NaN;

%% Plot projection

if isDes % make electron velocities 10^3 km/s
    vUnitStr= '(10^3 km/s)'; vUnitFactor = 1e-6;
else % ion velocities km/s
    vUnitStr= '(km/s)'; vUnitFactor = 1e-3;
end

pcolor(ax,vx_mesh*vUnitFactor,vy_mesh*vUnitFactor,log10(Fg')); 
shading(ax,'flat')

ax.XLabel.String = [vlabelx ' ' vUnitStr];
ax.YLabel.String = [vlabely ' ' vUnitStr];
ax.ZLabel.String = [vlabelz ' ' vUnitStr];

% set clim in whatever units
if have_clim; ax.CLim = clim; end

% set vlim in km
if have_vlim, ax.YLim = vlim*[-1 1]*vUnitFactor*1e3; pause(.01); ax.XLim = ax.YLim; end

hold(ax,'on')
axis(ax,'equal')

if doFlipX; ax.XDir = 'reverse'; end
if doFlipY; ax.YDir = 'reverse'; end

title(ax,dist.time(it).toUtc)

% last output
if nargout == 4
   plspec = []; plspec.F = Fg; plspec.vx = vx_mesh; plspec.vy = vy_mesh; 
end

end



function [Fg,vxMesh,vyMesh,dens,velmom] = project_sph_distr(F,v,phi,th,zphat,xphat,vg,phig,nMC,vzint)
%PROJECT_SPH_DISTR The function that performs the integration

% complete RH system
yphat = cross(zphat,xphat);

% diffs
dV = diff(v); dV = [dV(1),dV]; % quick and dirty
dPhi = abs(median(diff(phi))); % constant
dTh = abs(median(diff(th))); % constant

% primed diffs
dVg = diff(vg); dVg = [dVg(1),dVg]; % quick and dirty
dPhig = median(diff(phig)); % constant

% Number of projection bins
nVg = length(vg);
nAzg = length(phig);

% bin edges
phig_edges = [phig-dPhig/2,phig(end)+dPhig/2];
vg_edges = [vg(1)-dVg(1)/2,vg+dVg/2]; % quick and dirty

% convert to cartesian
[phiMesh,vMesh] = meshgrid(phig_edges+dPhig/2,vg); % Creates the mesh
[vxMesh,vyMesh] = pol2cart(phiMesh-pi/nAzg,vMesh);    % Converts to cartesian


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
% Area element (primed)
dAg = vg.*dVg*dPhig;

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
            vxp = sum([vx,vy,vz].*xphat,2);
            vyp = sum([vx,vy,vz].*yphat,2);
            vzp = dot([vx(1),vy(1),vz(1)],zphat);
            
            if vzp < vzint(1) || vzp > vzint(2); F(i,j,k) = 0; end
            
            % convert to polar coordinates (phip could become negative)
            [phip,vp] = cart2pol(vxp,vyp);
            
            % not so good but better than throwing away data? Should be
            % very rare anyway.
            vp(vp<vg_edges(1)*1.01) = vg_edges(1)*1.01;
            vp(vp>vg_edges(end)*.99) = vg_edges(end)*.99;
            
            % fix if negative
            phip(phip<0) = 2*pi+phip(phip<0);
            
            % get best indices. these two lines are the slowest in the
            % function so they should be optimized first!
            iAzg = interp1([phig(1)-dPhig,phig,phig(end)+dPhig],0:nAzg+1,phip,'nearest');
            iVg = interp1([vg(1)-dVg(1),vg,vg(end)+dVg(end)],0:nVg+1,vp,'nearest');
            
            % Loop through MC points and add value of instrument bin to the
            % appropriate projection bin
            for l = 1:nMC
                try
                    Fg(iAzg(l),iVg(l)) = Fg(iAzg(l),iVg(l))+F(i,j,k)*dtau(i,j,k)/dAg(iVg(l))/nMC;
                catch
                    disp(['Something went wrong. iAzg = ',num2str(iAzg(l)),',  iVg = ',num2str(iVg(l))])
                end
            end
        end
    end
end

% fix for interp shading, otherwise the last row can be whatever
Fg(end,:) = mean([Fg(end-1,:);Fg(1,:)]);

% Calculate density if requested
if nargout >= 4
    dAG = repmat(dAg,nAzg,1);
    dens = nansum(nansum(Fg(1:end-1,:).*dAG));
end
% Calculate velocity moment if requested
if nargout >= 5
    VG = repmat(vg',1,nVg);
    PHIG = repmat(phig,nAzg,1);
    [VXG,VYG] = pol2cart(PHIG,VG);
    
    velmom = [0,0];
    for l = 1:nVg
        for m = 1:nAzg
            velmom = velmom+[VXG(l,m),VYG(l,m)]*Fg(m,l)*dAg(l);
        end
    end
    velmom = velmom/dens;
end

end

