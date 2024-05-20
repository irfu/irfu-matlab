function [ax,hcb,C] = plot_skymap(varargin)
% MMS.PLOT_SKYMAP Plots skymap.
%
%  [ax,hcb] = MMS.PLOT_SKYMAP(PDist,'Opt1',OptVal1,...);
%     ax - handle to axes
%     hcb - handle to colorbar
%     C - skymap data
%  MMS.PLOT_SKYMAP(AX,...); - plot in axes AX.
%
%  Options:
%    'tint' - plot data for time interval if tint.length = 2
%             or closest time if tint.length = 1
%    'energy' - plot energy closest to given energy
%    'energylevel' - plot given energylevel
%    'vectors' - Nx2 cell array with 1x3 vector in first column and
%                textlabel in second column, eg.
%                vectors = {Bhat,'B';Ehat,'E'}
%    'vectorlabelcolor' - 'k', 'r', [1 1 1];
%    'flat' - plot a flat skymap with polar angle 0 at the top (ie. not a sphere)
%    'log' - plot log10 scale
%    'energytable' - energytable from v1 data
%    'phi' - phi data froom V1 data, vary from time to time, ax.XLim = [-10 370];
%    'phib' & 'polarb' [no theta] - solid angle boundary, only for 'flat';
%    'normal' - plot 3D circle with normal; 2020-03-24, wy
%    'avg_reduce_integral' - 'avg': distribution function average; [default]
%                          - 'reduce': f * delta V;
%                          - 'integral' f * Velocity Volumn; 2020-03-24;

[ax,args,nargs] = axescheck(varargin{:});

dist = args{1}; args = args(2:end);
if isempty(dist); irf.log('warning','Empty input.'); return; end

doSmooth = 0;
plotLog = 0;
fString = ['(' dist.units ')'];
plotSphere = 1;
plot3DCircle = 0;               % keyword for 'normal';
plotb = 0;
flag_energy = 0;
have_vectors = 0;
avg_reduce_integral = 'avg';
vectorlabelcolor = [1 1 1];
tId = 1:dist.length;
eId = 1:32;
tint = dist.time([1 end]);
[~,xi] = hist([log10(10),log10(30e3)],32); energyTable = 10.^xi;
% Set up spherical coordinate system.
r = 1; % radius of sphere
phi_edges = linspace(0,2*pi,size(dist.data,3)+1);  % azimuthal angle bin edges, default
theta_edges = linspace(0,pi,size(dist.data,4)+1); % polar angle bin edges, default
doPitchangles = 0;
doDouble = 0;

if nargs > 1, have_options = 1; else have_options = 0; end

while have_options
  l = 1;
  switch(lower(args{1}))
    case 'energytables'
      l = 2;
      energyTable = args{2};
    case 'energy'
      l = 2;
      energy = args{2};
      flag_energy = 1;
      eId = find(abs(energyTable-energy)==min(abs(energyTable-energy)));
    case 'energylevel'
      l = 2;
      eId = args{2};
    case 'phi_edges'            % phi edges from V1/V2 data
      l = 2;
      phi_edges = args{2};
    case 'phib'            % phi boundary for picking partial distribution
      l = 2;
      phib = args{2};
      plotb = 1;            % flag --> 1
    case 'polarb'            % theta boundary for picking partial distribution
      l = 2;
      polarb = args{2};
    case 'tint'
      l = 2;
      tint = args{2};
      if tint.length == 1
        tId = find(abs(dist.time-tint)==min(abs(dist.time-tint)));
      else
        [tId,~] = dist.time.tlim(tint);
        if isempty(tId); irf.log('warning','No data for given time interval.'); return; end
      end
    case 'vectors'
      l = 2;
      vectors = args{2};
      have_vectors = 1;
    case 'vectorlabelcolor'
      l = 2;
      vectorlabelcolor = args{2};
    case {'pitchangle','pitchangles'}
      l = 3;
      B_for_pitchangles = args{2};
      pitchangle_levels = args{3};
      doPitchangles = 1;
    case 'flat'
      plotSphere = 0;
    case 'sphere'
      plotSphere = 1;
    case 'double'
      doDouble = 1;
    case 'normal'       % normal vector;
      l = 2;
      normal = args{2};
      plot3DCircle = 1;
    case 'avg_reduce_integral'
      l = 2;
      avg_reduce_integral = args{2};
    case {'log'}
      if plotLog ~= 1
        plotLog = 1;
        fString = ['log_{10}' fString];
      end
    case 'smooth'
      l = 2;
      nSmooth = args{2};
      doSmooth = 1;
  end
  args = args(l+1:end);
  if isempty(args), break, end
end

if strcmp(dist.type, 'skymap')
  energyTable_all = dist.depend{1};
  energyTable = irf.nanmean(energyTable_all(tId, :), 1);
  if flag_energy
    eId = find(abs(energyTable-energy)==min(abs(energyTable-energy)));
  else
    eId = 1:numel(energyTable);
  end
end

[PHI,THETA] = meshgrid(phi_edges,theta_edges);
X = -r*sin(THETA).*cos(PHI); % '-' because the data shows which direction the particles were coming from
Y = -r*sin(THETA).*sin(PHI);
Z = -r*cos(THETA);
% dist.data has dimensions nT x nE x nAz x nPol
Units = irf_units;
dist_size = size(dist.data);
switch(avg_reduce_integral)
  case 'avg'
    C = squeeze(mean(mean(dist.data(tId,eId,:,:),2,'omitnan'),1,'omitnan'))';
  case 'reduce'
    dist_tmp = dist.data(tId,eId,:,:);
    delta_energy_minus = dist.ancillary.delta_energy_minus(tId, :);
    delta_energy_plus = dist.ancillary.delta_energy_plus(tId, :);
    tmp = dist.energy;
    energy_plus = tmp(tId, :) + delta_energy_plus;
    energy_minus = tmp(tId, :) - delta_energy_minus;
    if strcmp(dist.species, 'electrons')
      vv_plus = sqrt(energy_plus * Units.e * 2/ Units.me);
      vv_minus = sqrt(energy_minus * Units.e * 2/ Units.me);
      delta_vv = vv_plus - vv_minus;
    elseif strcmp(dist.species, 'ions')
      vv_plus = sqrt(energy_plus * Units.e * 2/ Units.mp);
      vv_minus = sqrt(energy_minus * Units.e * 2/ Units.mp);
      delta_vv = vv_plus - vv_minus;
    end
    delta_vv = repmat(delta_vv, 1, 1, dist_size(3), dist_size(4));
    delta_vv = delta_vv(:, eId, :, :);
    C = dist_tmp .* delta_vv * 1e12;            % [s^2/m^5]
    C = squeeze(mean(sum(C, 2, 'omitnan'),1,'omitnan'))';
    fString = '(s^2/m^5)';
  case 'integral'
    dist_tmp = dist.data(tId,eId,:,:);
    delta_energy_minus = dist.ancillary.delta_energy_minus(tId, :);
    delta_energy_plus = dist.ancillary.delta_energy_plus(tId, :);
    tmp = dist.energy;
    energy_plus = tmp(tId, :) + delta_energy_plus;
    energy_minus = tmp(tId, :) - delta_energy_minus;
    if strcmp(dist.species, 'electrons')
      vv_plus = sqrt(energy_plus * Units.e * 2/ Units.me);
      vv_minus = sqrt(energy_minus * Units.e * 2/ Units.me);
      delta_vv = vv_plus - vv_minus;
      vv2 = dist.ancillary.energy1 * 2 * Units.e / Units.me;      % bug here: assume energy0 = energy1
    elseif strcmp(dist.species, 'ions')
      vv_plus = sqrt(energy_plus * Units.e * 2/ Units.mp);
      vv_minus = sqrt(energy_minus * Units.e * 2/ Units.mp);
      delta_vv = vv_plus - vv_minus;
      vv2 = dist.ancillary.energy1 * 2 * Units.e / Units.mp;
    end
    delta_vv = repmat(delta_vv, 1, 1, dist_size(3), dist_size(4));
    delta_vv = delta_vv(:, eId, :, :);
    vv2 = repmat(vv2, length(tId), 1, dist_size(3), dist_size(4));
    sintheta = sind(dist.depend{3});
    sintheta = repmat(sintheta, length(tId), 1, dist_size(2), dist_size(3));
    sintheta = permute(sintheta,[1 3 4 2]);
    C = dist_tmp .* delta_vv .* vv2(:, eId, :, :) .* sintheta(:, eId, :, :) * 11.25/180 * pi * 11.25/180 * pi * 1e6;
    C = squeeze(mean(sum(C, 2, 'omitnan'), 1, 'omitnan'))';
    fString = '(cm^{-3})';
end

% Plot skymap
if isempty(ax), fig = figure; ax = axes; end
if doSmooth
  %smC = C;
  %smC()
  C = smooth2(double(C),nSmooth);
end
if plotLog, C = log10(C); end % units are whatever the input units were

if plotSphere
  hs = surf(ax,X,Y,Z,C);
  axis(ax,'square')
  axis(ax,'equal')
  ax.XLabel.String = 'X';
  ax.YLabel.String = 'Y';
  ax.ZLabel.String = 'Z';
  shading(ax,'flat');
else % plot flat map
  % change matrix so it corresponds to where the particles are going to, not coming from
  plotC = flipdim([C(:,17:32) C(:,1:16)],1);
  if doDouble
    plotC = [plotC,plotC];
    PHI_ = [PHI, 2*pi + PHI(:,2:end)]*180/pi;
    THETA_ = [THETA,THETA(:,2:end)]*180/pi;
    hs = surf(ax,PHI_,THETA_,THETA_*0,smooth2(double(plotC),1));
    %hs = surf(ax,PHI_,THETA_,THETA_*0,plotC);
  else
    hs = surf(ax,PHI*180/pi,THETA*180/pi,THETA*0,plotC);
  end
  ax.XLabel.String = 'Azimuthal angle (deg)';
  ax.YLabel.String = 'Polar angle (deg)';
  shading(ax,'flat');
  view(ax,[0 0 -1])
  ax.YLim = [0 180];
  ax.XLim = [0 360];
  if doDouble, ax.XLim = [0 720]; end
  ax.Box = 'on';
end
hcb = colorbar('peer',ax);
hcb.YLabel.String = fString;
titleString = {[irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'],['Energy = ' num2str(energyTable(eId),' %.0f') ,' eV']};
ax.Title.String = titleString;

% Plot vectors
hold(ax,'on');

while have_vectors
  if plotSphere
    vecHat = vectors{1,1}/norm(vectors{1,1});
    vecTxt = vectors{1,2};

    scale = 1.5;
    quiver3(ax,-scale*vecHat(1),-scale*vecHat(2),-scale*vecHat(3),vecHat(1),vecHat(2),vecHat(3),2*scale,'linewidth',2)
    scale = 1.7;
    axes(ax)
    text(double(scale*vecHat(1)),double(scale*vecHat(2)),double(scale*vecHat(3)),vecTxt,'fontsize',14)
    if plot3DCircle
      theta_circle=0:0.01:2*pi;
      center = [0, 0, 0];
      radius = 1;
      v = null(normal);
      points = repmat(center',1,size(theta_circle,2))+radius*(v(:,1)*cos(theta_circle)+v(:,2)*sin(theta_circle));
      plot3(ax, points(1,:), points(2,:), points(3,:), 'r-');
    end
  else % plot flat skymap
    vecHat = vectors{1,1}/norm(vectors{1,1});
    vecTxt = vectors{1,2};
    [azim,elev,r] = cart2sph(vecHat(1),vecHat(2),vecHat(3)); % tip of arrow
    if azim<0, azim = azim + 2*pi; end % from [-180 180] to [0 360]
    if elev<0; pol = pi/2 + abs(elev); else, pol = pi/2 - elev; end % from elevation to polar

    plot(ax,azim*180/pi,pol*180/pi,'o','linewidth',2,'markersize',12,'color',[1 0 0])
    plot(ax,azim*180/pi,pol*180/pi,'o','linewidth',0.5,'markersize',2,'color',[1 0 0],'markerfacecolor',[0 0 0])
    axes(ax); text(double(azim*180/pi),double(pol*180/pi),['   ' vecTxt],'fontsize',14,'HorizontalAlignment','left', 'color', vectorlabelcolor)

    [azim,elev,r] = cart2sph(-vecHat(1),-vecHat(2),-vecHat(3)); % back of arrow
    if azim<0, azim = azim + 2*pi; end % from [-180 180] to [0 360]
    if elev<0; pol = pi/2 + abs(elev); else, pol = pi/2 - elev; end % from elevation to polar

    plot3(ax,azim*180/pi,pol*180/pi,0,'o','linewidth',2,'markersize',12,'color',[1 0 0])
    plot3(ax,azim*180/pi,pol*180/pi,0,'x','linewidth',2,'markersize',12,'color',[1 0 0])
    axes(ax); text(double(azim*180/pi),double(pol*180/pi),['   ' vecTxt],'fontsize',14,'HorizontalAlignment','left', 'color', vectorlabelcolor)
  end
  vectors = vectors(2:end,:);
  if isempty(vectors), break, end
end

if doPitchangles
  if isa(B_for_pitchangles,'TSeries')
    Btmp = B_for_pitchangles.tlim(dist.time([1 dist.nt]) + 0.5*0.03*[-1 1]).data;
    btmp = mean(Btmp,1);
    btmp = btmp/sqrt(sum(btmp.^2));
  else
    btmp = B_for_pitchangles/sqrt(sum(B_for_pitchangles.^2));
  end
  dotPA = X*btmp(1) + Y*btmp(2) + Z*btmp(3);
  PA = acosd(dotPA);
  if numel(pitchangle_levels) == 1
    pitchangle_levels = pitchangle_levels*[1 1];
  end

  hold(ax,'on')
  climtmp = ax.CLim;
  if plotSphere

  else
    if doDouble
      plotPA = [PA,PA(:,2:end)];
      PHI_ = [PHI, 2*pi + PHI(:,2:end)];
      THETA_ = [THETA,THETA(:,2:end)];
      contour(ax,PHI_*180/pi,THETA_*180/pi,plotPA,pitchangle_levels,'k','linewidth',1)
    else
      contour(ax,PHI*180/pi,THETA*180/pi,PA,pitchangle_levels,'k','linewidth',1)
    end
  end
  ax.CLim = climtmp;
  hold(ax,'off')
end

if plotb
  if not(plotSphere)
    phib = [phib, phib(1)];
    polarb = [polarb, polarb(1)];
    plot(ax, phib, polarb, 'linewidth',2, 'color',[1 1 1])
  end
end
hold(ax,'off');

end