function [hax,hcb, FF] = plot_projection(varargin)
%MMS.PLOT_PROJECTION Plots projection on a specified plane.
%
% [ax,hcb] = MMS.PLOT_PROJECTION(dist,'Opt1',OptVal1,...) For PDist format
% [ax,hcb] = MMS.PLOT_PROJECTION(dist,phi,theta,stepTable,energy0,energy1,'Opt1',OptVal1,...)
% MMS.PLOT_PROJECTION(AX,...) - plot in axes AX.
%
% Input:
%       pdist - Electron or ion distribution as a TSeries
%       phi - Phi angles as a TSeries
%       theta - theta angles as a vector or structure
%       stepTable - Energy step table associated with the distribution (TSeries)
%       energy0 - energy table 0 as a structure or vector
%       energy1 - energy table 1 as a structure or vector
%
% Options:
%   'tint' - plot data for time interval if tint.length = 2
%             or closest time if tint.length = 1. For tint includes two or
%             more distributions the energies are rebinned into 64
%             channels.
%   'vectors' - Nx2 cell array with 1x3 vector in first column and
%                textlabel in second column
%   'xyz' - 3x3 matrix with [x;y;z]. z is normal to the plotted plane and
%           x and y are made orthogonal to z and each other if they are
%           not already. If you want to plot different planes you have to
%           rotate this matrix -> [y;z;x] -> [z;x;y]
%   'clim' - [cmin cmax], colorbar limits in logscale
%   'vlim' - vmax in km/s, zoom in to XLim = YLim = vlim*[-1 1]
%   'elevationlim' - in degrees [0 90], limiting elevation angle
%                    above/below projection plane to include in projection
%   'scpot' - scpot (TSeries), Correct velocities for spacecraft potential. For a single
%   value of tint the closest potential is used. For an interval the
%   spacecraft potential is average over that interval.
%
%   Notes: If phi, theta, stepTable, energy0, energy1 are not specified
%          then assumed values are used.
%          If the final distribution to be plotted is averaged over two or
%          more distributions the final distribution is rebinned into 64
%          energy channels.
%
%   Example:
%     mms.plot_projection(disDist1,'tint',tint+0.1*[-1 1],'xyz',[-1 0 0; 0 -1 0;hatB0],'vlim',1000);
%     mms.plot_projection(desDist1,'tint',tint,'xyz',[-1 0 0; 0 -1 0;hatB0],'vlim',40000);

[ax,args,nargs] = axescheck(varargin{:});
hax = ax;

irf.log('warning','Consider using mms.plot_int_projection which integrates the distribution and has an improved binning method.');
irf.log('warning','Please verify that you think the projection is done properly!');

anglespassed = 0;
lengthE = 32;
notint = 1;

dist = args{1};
% Check if it's electron or ions
if strfind(dist.name,'des')
  isDes = 1;
elseif strfind(dist.name,'dis')
  isDes = 0;
else
  irf.log('warning','Can''t recognize if input is electron or ions. Assuming it''s electrons.');
  isDes = 1;
end

distunits = 's^3cm^{-6}'; % These used are assumed if data is not PDist format

if isempty(dist); irf.log('warning','Empty input.'); return; end

% Check inputs for angles and energies
if nargs>1 && isa(args{2},'TSeries')
  irf.log('notice','Angles and energies passed.')
  phi = args{2};
  theta = args{3};
  stepTable = args{4};
  energy0 = args{5};
  energy1 = args{6};
  args = args(7:end);
  if isstruct(theta), theta = theta.data; end
  if isstruct(energy0), energy0 = energy0.data; end
  if isstruct(energy1), energy1 = energy1.data; end
  polar = theta*pi/180;
  dE = median(diff(log10(energy0)))/2;
  energy0Edges = 10.^(log10(energy0)-dE);
  energy0Edges = [energy0Edges 10.^(log10(energy0(end))+dE)];
  energy1Edges = 10.^(log10(energy1)-dE);
  energy1Edges = [energy1Edges 10.^(log10(energy1(end))+dE)];
  anglespassed = 1;
else
  if isa(dist,'PDist')
    irf.log('notice','Data is PDist format.')
    %Check if data is skymap
    if ~strcmp(dist.type,'skymap')
      irf.log('critical','PDist must be skymap format.');
      return;
    end
    tmpPhi = dist.depend{1,2};
    if size(tmpPhi,1) == 1 % only one value for one time or for all times
      phi = TSeries(dist.time,repmat(dist.depend{1,2},dist.length,1));
    else
      phi = TSeries(dist.time,dist.depend{1,2});
    end
    theta = dist.depend{1,3};
    polar = theta*pi/180;
    stepTable = TSeries(dist.time,dist.ancillary.esteptable);
    if and(isfield(dist.ancillary, 'energy0'), isfield(dist.ancillary, 'energy1'))
      energy0 = dist.ancillary.energy0;
      energy1 = dist.ancillary.energy1;
    else
      if isfield(dist.ancillary, 'energy')
        energy0 = dist.ancillary.energy(1, :);
        energy1 = dist.ancillary.energy(2, :);
      else
        irf.log('warning', 'no data for energy0 & energy1.');
      end
    end
    distunits = dist.units;
    dE = median(diff(log10(energy0)))/2;
    energy0Edges = 10.^(log10(energy0)-dE);
    energy0Edges = [energy0Edges 10.^(log10(energy0(end))+dE)];
    energy1Edges = 10.^(log10(energy1)-dE);
    energy1Edges = [energy1Edges 10.^(log10(energy1(end))+dE)];
    anglespassed = 1;
  else
    irf.log('warning','Angles and energies not passed. Assumed values used.')
    energyEdges = 10.^linspace(log10(10),log10(30e3),33);
    % Set up spherical coordinate system.
    [~,polar] = hist([0 pi],16);
    [~,azimuthal] = hist([0,2*pi],32);
  end
  args = args(2:end);
end

doFlipX = 0;
doFlipY = 0;
have_vlabels = 0;
have_vectors = 0;
have_clim = 0;
have_vlim = 0;
limElevation = 20;
correctForBinSize = 0;
includescpot = 0;
changed_xyz = [0 0 0];

x = [1 0 0]; y = [0 1 0]; z = [0 0 1]; % default vectors

have_options = 0;
if nargs > 1, have_options = 1; end
while have_options
  l = 1;
  switch(lower(args{1}))
    case 'tint'
      l = 2;
      notint = 0;
      tint = args{2};
      if tint.length == 1
        [~,tId] = min(abs(dist.time-tint));
        dist = TSeries(dist.time(tId),dist.data(tId,:,:,:));
        if anglespassed
          stepTable = TSeries(stepTable.time(tId),stepTable.data(tId));
          azimuthal = squeeze(phi.data(tId,:))*pi/180;
          if stepTable.data
            energyEdges = energy1Edges;
          else
            energyEdges = energy0Edges;
          end
        end
      else
        dist = dist.tlim(tint);
        if (length(dist.time)<1); irf.log('warning','No data for given time interval.'); return; end
        if anglespassed
          stepTable = stepTable.tlim(tint);
          phi = phi.tlim(tint);
          if (length(dist.time) > 1)
            irf.log('notice','Rebinning distribution.')
            [dist,phi,energy] = mms.psd_rebin(dist,phi,energy0,energy1,stepTable);
            lengthE = 64;
            azimuthal = phi.data*pi/180;
            energyEdges = 10.^(log10(energy)-dE/2);
            energyEdges = [energyEdges 10.^(log10(energy(end))+dE/2)];
          else
            azimuthal = phi.data*pi/180;
            if stepTable.data
              energyEdges = energy1Edges;
            else
              energyEdges = energy0Edges;
            end
          end
        end
      end
    case 'vectors'
      l = 2;
      vectors = args{2};
      have_vectors = 1;
    case 'xyz'
      l = 2;
      coord_sys = args{2};

      if ~all(size(coord_sys) == [3 3]) % to avoid errors associated with coordinate system: check for size of xyz, give error if not 3x3
        error('Please verify your coordinate system ''xyz'', it has to be a 3x3 matrix: [x;y;z], where x,y,z are 1x3 vectors.')
      end

      x = coord_sys(1,:)/norm(coord_sys(1,:));
      y = coord_sys(2,:)/norm(coord_sys(2,:));
      z = coord_sys(3,:)/norm(coord_sys(3,:));
      z = cross(x,y); z = z/norm(z);
      y = cross(z,x); y = y/norm(y);

      if abs(acosd(x*(coord_sys(1,:)/norm(coord_sys(1,:)))'))>1 % x is never changed
        irf.log('warning',['In making ''xyz'' a right handed orthogonal coordinate system, x (in-plane 1) was changed from [' num2str(coord_sys(1,:)/norm(coord_sys(1,:)),'% .2f') '] to [' num2str(x,'%.2f') ']. Please verify that this is according to your intentions.']);
        changed_xyz(1) = 1;
      end
      if abs(acosd(y*(coord_sys(2,:)/norm(coord_sys(2,:)))'))>1
        irf.log('warning',['In making ''xyz'' a right handed orthogonal coordinate system, y (in-plane 2) was changed from [' num2str(coord_sys(2,:)/norm(coord_sys(2,:)),'% .2f') '] to [' num2str(y,'%.2f') ']. Please verify that this is according to your intentions.']);
        changed_xyz(2) = 1;
      end
      if abs(acosd(z*(coord_sys(3,:)/norm(coord_sys(3,:)))'))>1
        irf.log('warning',['In making ''xyz'' a right handed orthogonal coordinate system, z (out-of-plane) was changed from [' num2str(coord_sys(3,:)/norm(coord_sys(3,:)),'% .2f') '] to [' num2str(z,'%.2f') ']. Please verify that this is according to your intentions.']);
        changed_xyz(3) = 1;
      end
    case 'clim'
      l = 2;
      clim = args{2};
      have_clim = 1;
    case 'vlim'
      l = 2;
      vlim = args{2};
      have_vlim = 1;
    case 'elevationlim'
      l = 2;
      limElevation = args{2};
    case 'vlabel'
      l = 2;
      have_vlabels = 1;
      vlabels = args{2};
    case 'usebincorrection'
      l = 2;
      correctForBinSize = args{2};
    case 'flipx'
      doFlipX = 1;
    case 'flipy'
      doFlipY = 1;
    case 'scpot'
      l = 2;
      scpot = args{2};
      if isa(scpot,'TSeries')
        includescpot = 1;
        irf.log('notice','Spacecraft potential passed.')
      else
        includescpot = 0;
        irf.log('notice','scpot not recognized. Not using it.')
      end
  end
  args = args(l+1:end);
  if isempty(args), break, end
end

if have_vlabels
  if changed_xyz(1)
    vlabelx = ['v_{x=[' num2str(x,'% .2f') ']}'];
    irf.log('warning',['x was changed, overriding given vlabels.']);
  else
    vlabelx = vlabels{1};
  end
  if changed_xyz(2)
    vlabely = ['v_{y=[' num2str(y,'% .2f') ']}'];
    irf.log('warning',['y was changed, overriding given vlabels.']);
  else
    vlabely = vlabels{2};
  end
  if changed_xyz(3)
    vlabelz = ['v_{z=[' num2str(z,'% .2f') ']}'];
    irf.log('warning',['z was changed, overriding given vlabels.']);
  else
    vlabelz = vlabels{3};
  end
else
  vlabelx = ['v_{x=[' num2str(x,'% .2f') ']}'];
  vlabely = ['v_{y=[' num2str(y,'% .2f') ']}'];
  vlabelz = ['v_{z=[' num2str(z,'% .2f') ']}'];
end

if (notint && anglespassed)
  if (length(dist.time) > 1)
    [dist,phi,energy] = mms.psd_rebin(dist,phi,energy0,energy1,stepTable);
    irf.log('notice','Recompling distribution into 64 energy channels.')
    lengthE = 64;
    azimuthal = phi.data*pi/180;
    energyEdges = 10.^(log10(energy)-dE/2);
    energyEdges = [energyEdges 10.^(log10(energy(end))+dE/2)];
  else
    azimuthal = phi.data*pi/180;
    if stepTable.data
      energyEdges = energy1Edges;
    else
      energyEdges = energy0Edges;
    end
  end
end

% Construct polar and azimuthal angle matrices
polar = ones(length(dist.time),1)*polar;
if (anglespassed==0)
  azimuthal = ones(length(dist.time),1)*azimuthal;
end
FF = zeros(length(dist.time),32,lengthE); % azimuthal, energy
edgesAz = linspace(0,2*pi,33);

for ii = 1:length(dist.time)
  [POL,AZ] = meshgrid(polar(ii,:),azimuthal(ii,:));
  X = -sin(POL).*cos(AZ); % '-' because the data shows which direction the particles were coming from
  Y = -sin(POL).*sin(AZ);
  Z = -cos(POL);

  % Transform into different coordinate system
  xX = reshape(X,size(X,1)*size(X,2),1);
  yY = reshape(Y,size(Y,1)*size(Y,2),1);
  zZ = reshape(Z,size(Z,1)*size(Z,2),1);

  newTmpX = [xX yY zZ]*x';
  newTmpY = [xX yY zZ]*y';
  newTmpZ = [xX yY zZ]*z';

  newX = reshape(newTmpX,size(X,1),size(X,2));
  newY = reshape(newTmpY,size(X,1),size(X,2));
  newZ = reshape(newTmpZ,size(X,1),size(X,2));

  elevationAngle = atan(newZ./sqrt(newX.^2+newY.^2));
  planeAz = atan2(newY,newX)+pi;

  % gets velocity in direction normal to 'z'-axis
  geoFactorElev = cos(elevationAngle);
  % geoFactorBinSize - detector bins in 'equator' plane are bigger and get a larger weight.
  % I think this is not good for the implementation in this function
  if correctForBinSize, geoFactorBinSize = sin(POL);
  else, geoFactorBinSize = 1; end

  % New coordinate system
  [nAz,binAz] = histc(planeAz,edgesAz);

  % Collect data in new distribution function FF
  for iE = 1:lengthE
    for iAz = 1:32
      % dist.data has dimensions nT x nE x nAz x nPol
      C = squeeze(dist.data(ii,iE,:,:));
      %C = squeeze(nansum(dist.data(tId,iE,:,:),1));
      C = C.*geoFactorElev.*geoFactorBinSize;
      C(abs(elevationAngle)>limElevation*pi/180) = NaN;
      % surf(ax,newX,newY,newZ,C);xlabel(ax,'x');ylabel(ax,'y')
      C(planeAz<(edgesAz(iAz)-0.1*pi/180)) = NaN;             % use 0.1 deg to fix Az angle edges bug
      C(planeAz>(edgesAz(iAz+1)+0.1*pi/180)) = NaN;           % use 0.1 deg to fix Az angle edges bug

      FF(ii,iAz,iE) = irf.nanmean(irf.nanmean(C));
    end
  end
end

if isscalar(dist.time)
  FF = squeeze(FF);
  tint = dist.time;
else
  FF = squeeze(irf.nanmean(FF,1));
  tint = irf.tint(dist.time.start.utc,dist.time.stop.utc);
end

% Plot projection plane
if isempty(ax), fig = figure; ax = axes; axis(ax,'square'); end

% Calculate spacecraft potential for distribution
if includescpot
  if isscalar(tint)
    [~,tId] = min(abs(scpot.time-tint));
    scpot = scpot.data(tId);
  else
    scpot = scpot  .tlim(tint);
    scpot = irf.nanmean(scpot.data);
  end
else
  scpot = 0;
end


Units = irf_units; % Use IAU and CODATA values for fundamental constants.
if isDes
  m = Units.me;
else
  m = Units.mp;
  scpot = -scpot;
end

speedTable = real(sqrt((energyEdges-scpot)*Units.e*2/m)*1e-3); % km/s

rE = speedTable;
%rE = energyTable;
plX = rE'*cos(edgesAz+pi);
plY = rE'*sin(edgesAz+pi);

FF(FF <= 0) = NaN; % set to white the zero points

if isDes % make electron velocities 10^3 km/s
  hs = surf(ax,plX*1e-3,plY*1e-3,plY*0,log10(FF'));
  vUnitStr= '(10^3 km/s)';
  if have_vlim, ax.YLim = vlim*[-1 1]*1e-3; ax.XLim = ax.YLim; end
  hcb = colorbar('peer',ax);
  hcb.YLabel.String = ['log_{10} f_e (' distunits ')'];
else % ion velocities km/s
  hs = surf(ax,plX,plY,plY*0,log10(FF'));
  vUnitStr= '(km/s)';
  if have_vlim, ax.YLim = vlim*[-1 1]; ax.XLim = ax.YLim; end
  hcb = colorbar('peer',ax);
  hcb.YLabel.String = ['log_{10} f_i (' distunits ')'];
end

ax.XLabel.String = [vlabelx ' ' vUnitStr];
ax.YLabel.String = [vlabely ' ' vUnitStr];
ax.ZLabel.String = [vlabelz ' ' vUnitStr];

axis(ax,'square')
axes(ax)
view([0 0 1])

if doFlipX
  cbPos = hcb.Position;
  axPos = ax.Position;
  view(ax,180,90)
  hcb.Position = cbPos;
  ax.Position = axPos;
end
ax.Box = 'on';

titleString = {[irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
ax.Title.String = titleString;

shading(ax,'flat');
if have_clim, ax.CLim = clim; end

% Plot vectors
hold(ax,'on');

while have_vectors
  if isempty(vectors), break, end
  vecHat = vectors{1,1}/norm(vectors{1,1});
  vecHat = vecHat;
  vecTxt = vectors{1,2};
  scale = ax.YLim(2)/3;
  scaleQuiver = 1.5*scale;
  quiver3(ax,0,0,0,(vecHat*x'),(vecHat*y'),abs(vecHat*z'),scaleQuiver,'linewidth',1,'linestyle','-','color','k')
  scaleTxt = 1.9*scale;
  text(double(scaleTxt*vecHat*x'),double(scaleTxt*vecHat*y'),1,vecTxt,'fontsize',10)
  vectors = vectors(2:end,:);
  axis(ax,'square')
end
hold(ax,'off');

end
