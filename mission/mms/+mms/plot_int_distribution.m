function [ax,pst] = plot_int_distribution(varargin)
% NOTICE:
% This function has been replaced by PDist.reduce and is no longer
% being updated.
%
%MMS.PLOT_INT_DISTRIBUTION Plots integrated FPI data along a specified vector.
%
%   Does not seem optimal for cuts in the distribution function. Works fine
%   for ion spectrograms.
%
%   MMS.PLOT_INT_DISTRIBUTION(dist,'Opt1',OptVal1,...) plots projection of
%   particle data in the PDist object dist. The plot can be a time series
%   spectrogram or a time-averaged line plot.
%
%   MMS.PLOT_INT_DISTRIBUTION(AX,...) plots in axes AX instead of current
%   axes.
%
%   [hsf,pst] = MMS.PLOT_INT_DISTRIBUTION(...) returns surface/line and
%   a structure pst with data used for plotting.
%   pst contains for a line:
%       F       - particle flux, for PSD, typically the unit is [s m^-4],
%               which 10^-3*[s km^-1 cm^-3]
%       v       - v grid in [m/s]
%       dens    - number density, time-averaged, derived from the
%               integrated data in [m^-3].
%       vel     - 1D velocity, time-averaged, derived from the integrated
%               data in [m/s].
%   plspec contains for a spectrogram:
%       p       - particle flux, for PSD, typically the unit is [s m^-4]
%       f       - v grid in [m/s]
%       t       - time vector
%       dens    - number density, as a vector, derived from the integrated
%               data in [m^-3]
%       vel     - 1D velocity, as a vector, derived from the integrated
%               data in [m/s]
%   To recreate plot elsewhere use:
%       Line:
%           plot(plspec.v,plspec.F,'-o');set(gca,'YScale','Log')
%       Time series:
%           irf_spectrogram(plspec)
%   Options:
%   'tint'/'t'  - time interval used for plot. Data is averaged for a line
%               and all data in interval is shown for a spectrogram
%   'xvec'      - normal vector to plot data against, can be a TSeries
%               object or a row vector
%   'plotmode'  - can be 'spec'/'ts' or 'line', default is 'line'
%   'clim'      - [cmin cmax], colorbar limits in logscale or ylim for line
%   'vlim'      - vmax in [km/s], zoom in to XLim = YLim = vlim*[-1 1]
%   'vint'      - set limits on the from-line velocity to get cut-like
%               distribution
%   'nMC'       - number of Monte Carlo iterations used for integration,
%               default is 100
%   'scpot'     - TSeries object with spacecraft potential in [V]
%   'weight'-   how the number of MC iterations per bin is weighted, can be
%               'none' (default), 'lin' or 'log'
%   'vlabel'    - string for axis label corresponding to v
%   'flipx'     - boolean value where 1 flips the x axis
%   'vg'        - array with center values for the projection velocity
%               grid in [km/s], determined by instrument if omitted
%
%   Examples:
%       % Example 1: Flat top distribution at shock
%       t1 = irf.time_array('2017-10-24T20:02:40',0);
%       t2 = irf.time_array('2017-10-24T20:03:15',0);
%       db_info = datastore('mms_db');
%       file  = [db_info.local_file_db_root,...
%           '/mms1/fpi/brst/l2/des-dist/2017/10/24/mms1_fpi_brst_l2_des-dist_20171024200103_v3.2.0.cdf'];
%       ePDist = mms.make_pdist(file);
%       figure;mms.plot_int_distribution(ePDist,'t',t1,'xvec',[-1,0,1]/sqrt(2),'nmc',200,'aint',20*[-1,1],'vg',linspace(2e3,20e3,32)); hold on;
%       mms.plot_int_distribution(ePDist,'t',t2,'xvec',[-1,0,1]/sqrt(2),'nmc',200,'aint',20*[-1,1],'vg',linspace(2e3,20e3,32));
%       set(gca,'XLim',[0,20]); set(gca,'XScale','log'); legend('Upstream','Downstream')
%
%       % Example 2: Ion phase space holes due to shock ripples
%       tint = irf.tint('2015-10-07T11:44:40/2015-10-07T11:44:45');
%       db_info = datastore('mms_db');
%       file  = [db_info.local_file_db_root,...
%           '/mms1/fpi/brst/l2/dis-dist/2015/10/07/mms1_fpi_brst_l2_dis-dist_20151007114414_v3.1.0.cdf'];
%       iPDist = mms.make_pdist(file);
%       mms.plot_int_distribution(iPDist,'t',tint,'plotmode','ts','xvec',[0.88 0.46 -0.11],'nmc',10,'vlim',800);
%
%
%   The function uses a Monte Carlo integration method. To read more about
%   it:
%
%   See also: IRF_INT_SPH_DIST MMS.PLOT_INT_PROJECTION

%   Written by: Andreas Johlander, andreasj@irfu.se
%
%   TODO:   Recalculate energies given spacecraft potential input
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
have_vlabels = 0;
have_vlim = 0;
have_clim = 0;
nMC = 100; % number of Monte Carlo iterations
vint = [-Inf,Inf];
aint = [-180,180];
plotMode = 'line';
vgInput = 0;
weight = 'none';
correctForScPot = 0;

have_options = nargs > 1;
while have_options
  switch(lower(args{1}))
    case {'t','tint'} % time
      tint = args{2};
    case 'xvec' % vector to plot data against
      xphat = args{2};
      if isnumeric(xphat) % only single vector input allowed
        xphat = xphat/norm(xphat);
        if size(xphat,2)==1; xphat = xphat'; end % make sure it's a row vector
      elseif isa(xphat,'TSeries') % time series
        xphat = xphat.resample(dist).norm.data;
      end

    case 'vlim'
      vlim = args{2};
      have_vlim = 1;
    case 'vlabel'
      have_vlabels = 1;
      vlabel = args{2};
    case 'clim'
      clim = args{2};
      have_clim = 1;
    case 'flipx'
      doFlipX = args{2};
    case 'nmc' % number of Monte Carlo iterations
      nMC = args{2};
    case 'vint' % limit on transverse velocity (like a cylinder) [km/s]
      vint = args{2};
    case 'aint'
      aint = args{2};
    case 'plotmode' % 'line' or 'spec'/'ts'
      plotMode = args{2};
    case 'vg' % define velocity grid
      vgInput = 1;
      vg = args{2}*1e3;
    case 'weight' % how data is weighted
      weight = args{2};
    case 'scpot'
      scPot = args{2};
      correctForScPot = 1;
  end
  args = args(3:end);
  if isempty(args), break, end
end

if ~have_vlabels
  if min(size(xphat)) == 1
    vlabel = ['v_{x=[' num2str(xphat,'% .2f') ']}'];
  else
    vlabel = 'v ';
  end
end




%% Get angles and velocities for spherical instrument grid, set projection
%  grid and perform projection

u = irf_units;

if isDes == 1; M = u.me; else; M = u.mp; end

% get all time indicies
if length(tint) == 1 % single time
  it = interp1(dist.time.epochUnix,1:length(dist.time),tint.epochUnix,'nearest');
else % time interval
  it1 = interp1(dist.time.epochUnix,1:length(dist.time),tint(1).epochUnix,'nearest');
  it2 = interp1(dist.time.epochUnix,1:length(dist.time),tint(2).epochUnix,'nearest');
  it = it1:it2;
end

% make sure xphat is a matrix
if min(size(xphat)) == 1
  xphat = repmat(xphat,length(dist),1);
end

% loop to get projection
for i = 1:length(it)
  if length(it)>1;disp([num2str(i),'/',num2str(length(it))]); end

  % 3d data matrix for time index it, [E,phi,th]
  F3d = double(squeeze(double(dist.data(it(i),:,:,:)))*1e12); % s^3/m^6

  emat = double(dist.energy); % in eV
  energy = emat(it(i),:);


  % if scpot is in input, correct for it
  if correctForScPot
    energy = energy-scPot.resample(dist).data(it(i));
    iEend = find(energy<0,1,'last');
    % set to zero if below sc potential, +1 to make sure no the bin is removed
    F3d(1:iEend+1,:,:) = 0;
    energy(energy<0) = 0; % not sure if needed
  end

  v = sqrt(2*energy*u.e/M); % m/s

  if length(v) ~= 32
    error('something went wrong')
  end

  % azimuthal angle
  phi = double(dist.depend{2}(it(i),:)); % in degrees
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


  % Set projection grid after the first distribution function
  if i == 1
    % bin centers
    if ~vgInput
      vg = [-fliplr(v),v];
    end
    % initiate projected f
    Fg = zeros(length(it),length(vg));
    dens = zeros(length(it),1);
    vel = zeros(length(it),1);
  end
  % perform projection
  tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat(it(i),:),'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight);
  Fg(i,:) = tmpst.F;
  dens(i) = tmpst.dens;
  vel(i) = tmpst.vel;
end


%% Plot distribution

if isDes % make electron velocities 10^3 km/s
  vUnitStr= '(10^3 km/s)'; vUnitFactor = 1e-6;
else % ion velocities km/s
  vUnitStr= '(km/s)'; vUnitFactor = 1e-3;
end

% Different plotting for line or time series
if strcmpi(plotMode,'line')
  % average
  Fg = nanmean(Fg,1);

  pst = [];
  pst.F = Fg;
  pst.v = vg;
  pst.dens = nanmean(dens,1);
  pst.vel = nanmean(vel,1);

  plot(ax,vg*vUnitFactor,Fg,'-o');
  ax.YScale = 'log';
  ax.XLabel.String = [vlabel ' ' vUnitStr];

  % set vlim in km
  if have_vlim, ax.XLim = vlim*[-1 1]*vUnitFactor*1e3; end

  % set clim
  if have_clim, ax.YLim = clim*[-1 1]; end

  if doFlipX; ax.XDir = 'reverse'; end


elseif strcmpi(plotMode,'ts') || strcmpi(plotMode,'spec')
  % do not average, instead initiate specrec for irf_spectrogram
  % put nans instead of 0s
  Fg(Fg==0) = NaN;
  pst = [];
  pst.t = dist.time(it).epochUnix;
  pst.f_label = 'velocity';
  pst.f = vg*vUnitFactor;
  pst.p = log10(Fg);

  irf_spectrogram(ax,pst);
  ax.YLabel.String = [vlabel ' ' vUnitStr];

  pst.dens = dens;
  pst.vel = vel;

  % set vlim in km
  if have_vlim, ax.YLim = vlim*[-1 1]*vUnitFactor*1e3; end

  % set clim
  if have_clim, ax.CLim = clim*[-1 1]; end

  if doFlipX; ax.YDir = 'reverse'; end
else
  error('Unknown plot mode')
end


% title is start time (=center time - dt/2) + time difference (=last center time + dt/2 - start time)
dt = median(diff(dist.time.epochUnix));
diffT = diff(dist.time([it(1),it(end)]).epochUnix+dt/2*[-1;1]);
title(ax,[dist.time(it(1)).toUtc,'+',num2str(diffT),'s'])

end
