function [hsf,pst] = plot_int_projection(varargin)
% NOTICE:
% This function has been replaced by PDist.reduce and is no longer
% being updated.
%
%MMS.PLOT_INT_PROJECTION Plots integrated FPI data on a specified plane.
%
%   MMS.PLOT_INT_PROJECTION(dist,'Opt1',OptVal1,...) plots projection of
%   particle data in the PDist object dist.
%
%   MMS.PLOT_INT_PROJECTION(AX,...) plots in axes AX instead of current
%   axes.
%
%   [hsf,plspec] = MMS.PLOT_INT_PROJECTION(...) returns surface
%   handle, and a structure plspec with data used for plotting.
%   plspec contains:
%       F   - the particle flux, for PSD, typically the unit is [s^2m^-5],
%           which is the same as [s^2cm^-3^km^-2]
%       vx      - vx grid in [m/s]
%       vy      - vy grid in [m/s]
%       dens    - number density, derived from the projected data in m^-3.
%       vel     - 2D velocity, derived from the projected data in m/s.
%
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
%               default is 500
%   'weight'-   how the number of MC iterations per bin is weighted, can be
%               'none' (default), 'lin' or 'log'
%   'vlabel'    - 1x3 cell array containing strings for axis labels
%               corresponding to x, y, and z
%   'flipx'/'flipy' - boolean value where 1 flips the x/y axis
%   'colorbar'  - boolean value where 1 adds a colorbar to plot
%   'vg'        - array with center values for the projection velocity
%               grid in [km/s], determined by instrument if omitted
%   'phig'      - array with center values for the projection azimuthal
%               grid in radians, determined by instrument if omitted
%
%   Examples:
%       tint = irf.tint('2015-10-16T13:07:02/2015-10-16T13:07:03');
%       t = irf.time_array('2015-10-16T13:07:02.226',0);
%       ePDist = mms.get_data('PDe_fpi_brst_l2',tint,4);
%       vg = linspace(-2.7e4,2.7e4,200);
%       mms.plot_int_projection(ePDist,'t',t,'vlim',1e4,'xyz',[0,1,0;-1,0,0; 0,0,1],'vzint',2000*[-1,1],'nmc',5e2,'vg',vg,'base','cart');
%
%
%   The function uses a Monte Carlo integration method. To read more about
%   it:
%
%   See also: IRF_INT_SPH_DIST

%   Written by: Andreas Johlander, andreasj@irfu.se
%
%   TODO:   Add time interval averaging
%           Recalculate energies given spacecraft potential input
%           Account for other quantities than psd, e.g. diff energy flux
%
%   Generally, the code deals with quantities in SI units until plotting.

%% Warning
irf.log('w','This function has been replaced by PDist.reduce')


%% Input

[ax,args,nargs] = axescheck(varargin{:});
if isempty(ax); ax = gca; end

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
nMC = 500; % number of Monte Carlo iterations
vzint = [-inf,inf];
showColorbar = 0;
vgInput = 0;
phigInput = 0;
weight = 'none';
base = 'pol';

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
    case 'colorbar' % 1 if show colorbar
      showColorbar = args{2};
    case 'vg' % define velocity grid
      vgInput = 1;
      vg = args{2}*1e3;
    case 'phig' % define velocity grid
      phigInput = 1;
      phig = args{2};
    case 'weight' % how data is weighted
      weight = args{2};
    case 'base'
      base = args{2};
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
it = interp1(dist.time.epochUnix,1:length(dist.time.epochUnix),t.epochUnix,'nearest');

% 3d data matrix for time index it
F3d = double(squeeze(double(dist.data(it,:,:,:)))*1e12); % s^3/m^6


emat = double(dist.energy); % in eV
energy = emat(it,:);
v = sqrt(2*energy*u.e/M); % m/s

if length(v) ~= 32
  error('something went wrong')
end

% azimuthal angle
phi = double(dist.depend{2}(it,:)); % in degrees
%phi = phi+180;
%phi(phi>360) = phi(phi>360)-360;
phi = phi-180;
phi = phi*pi/180; % in radians

if length(phi) ~= 32
  error('something went wrong')
end

% elevation angle
th = double(dist.depend{3}); % polar angle in degrees
th = th-90; % elevation angle in degrees
th = th*pi/180; % in radians

if length(th) ~= 16
  error('something went wrong')
end


%% Set projection grid
nAzg = 32;

% diffs
dPhig = 2*pi/nAzg;

% bin centers defined by user or set to same as instrument
if ~phigInput; phig = linspace(0,2*pi-dPhig,nAzg)+dPhig/2; end
if ~vgInput; vg = v; end % same as instrument if no input

%% perform projection
pst = irf_int_sph_dist(F3d,v,phi,th,vg,'z',zphat,'x',xphat,'phig',phig,'nMC',nMC,'vzint',vzint*1e3,'weight',weight,'base',base);

% put nans instead of 0s
pst.F(pst.F==0) = NaN;

%% Plot projection

if isDes % make electron velocities 10^3 km/s
  vUnitStr= '(10^3 km/s)'; vUnitFactor = 1e-6;
else % ion velocities km/s
  vUnitStr= '(km/s)'; vUnitFactor = 1e-3;
end

hsf = pcolor(ax,pst.vx*vUnitFactor,pst.vy*vUnitFactor,log10(pst.F'));
shading(ax,'flat')

ax.XLabel.String = [vlabelx ' ' vUnitStr];
ax.YLabel.String = [vlabely ' ' vUnitStr];
ax.ZLabel.String = [vlabelz ' ' vUnitStr];

% set clim in whatever units
if have_clim; ax.CLim = clim; end

% set vlim in km
if have_vlim, ax.YLim = vlim*[-1 1]*vUnitFactor*1e3; pause(.01); ax.XLim = ax.YLim; end

% show colorbar with appropriate units
if showColorbar
  fVarStr = 'log_1_0 F ';
  % chose unit string (add more if needed)
  if strcmp(dist.units,'s^3/cm^6'); fUnitStr = '[s^2/m^5]';
  else; fUnitStr = '';
  end
  hcb = colorbar(ax);
  ylabel(hcb,[fVarStr,fUnitStr])
end

hold(ax,'on')
axis(ax,'square')

if doFlipX; ax.XDir = 'reverse'; end
if doFlipY; ax.YDir = 'reverse'; end

title(ax,dist.time(it).toUtc)

end
