% Calculates the number of bowshock crossing based on the actual position
% of THOR and the bowshock location at the time.

%% Load THOR orbit
%datastore('spice','dir','/Users/Cecilia/calc/SPICE');
orbitKernels = 1;
resampleKernels = 0;
if orbitKernels
  units = irf_units;
  %rTHOR = thor_orbit('alt1a.bsp',2*3600);
  rTHOR = thor_orbit('new1a.bsp',2*3600);
  if resampleKernels
    newTime = rTHOR.time.start:120:rTHOR.time.stop; % 2 min intervals
    tmpR = rTHOR.resample(newTime);
    rTHOR = tmpR;
  end
else
  get_orbit
  rTHOR = irf.ts_vec_xyz(irf_time(t,'epoch>epochtt'),[x y x*0])
end

%% Download/load OMNI database data
% 3.3 years (duration of the orbit) of representative solar wind conditions
% should be chosen
loadDataFromFile = 1;
if loadDataFromFile
  load /Users/Cecilia/MATLAB/irfu-matlab/mission/thor/orbit_coverage/omni_data_20010101-20041231_bsnx_Bxyz_MS.mat
  % load /Users/Cecilia/MATLAB/irfu-matlab/mission/thor/orbit_coverage/omni_data_20090101-20121231_bsnx_Bxyz_MS.mat
  % load /Users/Cecilia/MATLAB/irfu-matlab/mission/thor/orbit_coverage/omni_data_20120101-20151231_bsnx_Bxyz_MS.mat
else % download data from omni database, change the years as is appropriate
  % Can only seem to download one years data at a time
  clear tintUTC
  tsub = 1;
  c_eval('tintUTC{tsub} = ''?-01-01T00:00:00/?-12-31T23:59:00''; tsub = tsub+1;',2001:2004);

  omni_orig = [];
  tic;
  for iy = 1:numel(tintUTC)
    tint = irf.tint(tintUTC{iy});
    tmp_omni = irf_get_data_omni(tint,'bsnx,Bx,By,Bz,Ms','omni_min');
    omni_orig = [omni_orig; tmp_omni];
    toc
  end

  % Clean up data
  omni = omni_orig;
  t0 = irf_time(omni_orig(:,1),'epoch>epochtt');

  if 0 % Removing all the points that dont have R0 data changes the total time
    omni(isnan(omni(:,2)),:)=[]; % remove all points that dont have R0 data
  end

  R0 = omni(:,2); % RE
  kmR0 = omni(:,2)*units.RE*1e-3; % km

  tsBSNX = irf.ts_scalar(irf_time(omni(:,1),'epoch>utc'),omni(:,2));
  tsBSNX = tsBSNX.resample(t0);
  tsBSNX.units = 'RE';
  tsBSNX.name = 'Bowshock nose distance, X';
  tsB = irf.ts_vec_xyz(irf_time(omni(:,1),'epoch>utc'),omni(:,3:5));
  tsB = tsB.resample(t0);
  tsB.units = 'nT';
  tsB.name = 'Solar wind magnetic field';
  tsM = irf.ts_scalar(irf_time(omni(:,1),'epoch>utc'),omni(:,6));
  tsM = tsM.resample(t0);
  tsM.units = '';
  tsM.name = 'Solar wind Mach number';
end

%% Adjust timelines of BSNX and THOR
units = irf_units;
t0 = irf_time(omni_orig(:,1),'epoch>epochtt');
%tR0 = irf_time(bsnx(:,1),'epoch>epochTT');
tShift=rTHOR.time.start-tsBSNX.time.start;
newTime = tsBSNX.time+tShift; % shift the time of bsnx to THOR's time
%newtR0 = newtR0.tlim(rTHOR.time([1 end]));

tsBSNXkm = irf.ts_scalar(newTime,tsBSNX.data*units.RE*1e-3); tsBSNXkm.units = 'km'; tsBSNXkm.name = 'Bowshock nose distance, X';
tsM = tsM.clone(newTime,tsM.data); tsM.units = ''; tsM.name = 'Solar wind Mach number';
tsB = irf.ts_vec_xyz(newTime,tsB.data); tsB.units = ''; tsB.name = 'Solar wind magnetic field';
tsBSNX = tsM.clone(newTime,tsBSNX.data); tsBSNX.units = 'RE'; tsBSNX.name = 'Bowshock nose distance, X';
%xBSN = irf.ts_scalar(newTime,tsBSNX); xBSN.units = 'km'; xBSN.name = 'Bowshock nose distance, X';

tsBSNXkm = tsBSNXkm.tlim(rTHOR.time([1 end]));
tsBSNX = tsBSNX.tlim(rTHOR.time([1 end]));
tsM = tsM.tlim(rTHOR.time([1 end]));
tsB = tsB.tlim(rTHOR.time([1 end]));
rTHOR = rTHOR.resample(tsBSNX); % upsample orbit times to bsnx's timeline, 1 min

%% Bowshock model
fr = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
fx = @(r,R0) 0.5*45.3/0.04-sqrt((0.5*45.3/0.04)^2+r.^2/0.04)+R0;

%% Check if THOR is inside bowshock or not
xTHOR = rTHOR.x.data/units.RE*1e3; % km->RE
yTHOR = rTHOR.y.data/units.RE*1e3;
zTHOR = rTHOR.z.data/units.RE*1e3;

r_THOR = sqrt(yTHOR.^2+zTHOR.^2);

rBS = fr(xTHOR,tsBSNX.data);
tsrBS = irf.ts_scalar(tsBSNX.time,rBS);

allInd = 1:rTHOR.length;
isInside = find(xTHOR<tsBSNX.data & abs(rBS)>abs(r_THOR));
isOutside = tocolumn(setdiff(allInd,isInside));

isInboundCrossing = isInside(find(diff(isInside)>1));
isOutboundCrossing = isOutside(find(diff(isOutside)>1));
isCrossing = sort([isOutboundCrossing; isInboundCrossing]);
nCrossings = numel(isCrossing);
%nCrossings = numel(isOutboundCrossing)+numel(isInboundCrossing);

nCrossingsPerYear = nCrossings/((rTHOR.time.stop-rTHOR.time.start)/60/60/24/365);
% Check time spent inside or outside the bowshock
dt = rTHOR.time(2)-rTHOR.time(1);
timeInside = dt*numel(isInside);
timeOutside = dt*numel(isOutside);

disp(sprintf('Dwell time, inside bowshock: %s days',num2str((timeInside)/60/60/24,'%g')))
disp(sprintf('Dwell time, outside bowshock: %s days',num2str((timeOutside)/60/60/24,'%g')))

%% Check quality of the bowshock crossing
Rout = 15;
QR = thor_QR(tsBSNX,Rout);
QV = thor_QV(tsBSNX,rTHOR);

% Using Bz=0 is mostly for bugcheck
tsB_Bz0 = tsB; tsB_Bz0.data(:,3) = 0;

[QBparBz0,AngleBparBz0,shockNormal] = thor_QB(tsBSNX,rTHOR,tsB_Bz0);
[QBpar,AngleBpar] = thor_QB(tsBSNX,rTHOR,tsB);
%[QBpar,AngleBpar,shockNormal] = thor_QB(tsBSNX(isCrossing),rTHOR(isCrossing),tsB(isCrossing));
[QBperpBz0,AngleBperpBz0] = thor_QB(tsBSNX,rTHOR,tsB_Bz0,'perp');
[QBperp,AngleBperp] = thor_QB(tsBSNX,rTHOR,tsB,'perp');


Qpar = QR*QBpar*QV;
Qperp = QR*QBperp*QV;

QparBz0 = QR*QBparBz0*QV;
QperpBz0 = QR*QBperpBz0*QV;

% Count crossings
nCrossBPerpBz0 = numel(find(QBperpBz0(isCrossing).data>0.8));
nCrossBParBz0 = numel(find(QBparBz0(isCrossing).data>0.8));
nCrossBPerp = numel(find(QBperp(isCrossing).data>0.8));
nCrossBPar = numel(find(QBpar(isCrossing).data>0.8));

nCrossQPerpBz0 = numel(find(QperpBz0(isCrossing).data>0.8));
nCrossQParBz0 = numel(find(QparBz0(isCrossing).data>0.8));
nCrossQPerp = numel(find(Qperp(isCrossing).data>0.8));
nCrossQPar = numel(find(Qpar(isCrossing).data>0.8));

edgesQ = 0:0.05:1;
nDistQpar = histc(Qpar(isCrossing).data,edgesQ);
nDistQBpar = histc(QBpar(isCrossing).data,edgesQ);
nDistQperp = histc(Qperp(isCrossing).data,edgesQ);
nDistQBperp = histc(QBperp(isCrossing).data,edgesQ);
nDistQV = histc(QV(isCrossing).data,edgesQ);

disp(sprintf('# Quasi-perp crossings, Bz=0, QB>0.8: %s',num2str(nCrossBPerpBz0,'%g')))
disp(sprintf('# Quasi-perp crossings, QB>0.8: %s',num2str(nCrossBPerp,'%g')))
disp(sprintf('# Quasi-perp crossings, Bz=0, Q>0.8: %s',num2str(nCrossQPerpBz0,'%g')))
disp(sprintf('# Quasi-perp crossings, Q>0.8: %s',num2str(nCrossQPerp,'%g')))

disp(sprintf('# Quasi-par crossings, Bz=0, QB>0.8: %s',num2str(nCrossBParBz0,'%g')))
disp(sprintf('# Quasi-par crossings, QB>0.8: %s',num2str(nCrossBPar,'%g')))
disp(sprintf('# Quasi-par crossings, Bz=0, Q>0.8: %s',num2str(nCrossQParBz0,'%g')))
disp(sprintf('# Quasi-par crossings, Q>0.8: %s',num2str(nCrossQPar,'%g')))

edgesAngle = [0:10:180];
centerAngle = edgesAngle(1:end)+edgesAngle(2)-edgesAngle(1);
nAngles = histc(AngleBpar.data(isCrossing,:),edgesAngle');
nAnglesBz0 = histc(AngleBparBz0.data(isCrossing,:),edgesAngle');

%% Figures: Set some figure parameters
markerStyle = 'o';
xlim = [-20 40];
ylim = [-30 30];
xlim = [-50 50];
ylim = [-50 50];
xtick = [-50:10:50];
ytick = [-50:10:50];
fontsize = 14;
plotB = 1;
plotShockNormal = 1;
cmap = colormap('jet');
mirrorcmap = [cmap; cmap(end:-1:1,:)];
colorOrbit = [0 0 0]+0.8;
plStep = 10;
doPrint = 0;
t1utc = t0(1).utc; t2utc = t0(end).utc;
time_str = sprintf('%s_%s',t1utc(1:10),t2utc(1:10));
printString = 'cn.print([time_str ''_box_'' num2str(diff(xlim)) ''x'' num2str(diff(ylim)) ''_'' figure_name],''thor'')';
scrsz = get(groot,'ScreenSize');

%% Figure: Shock normal angle in 3D
figure_name = 'thetaB_3D';
%figure_position = scrsz; figure_position(3) = figure_position(3)*0.5;
%figure('Position',figure_position)

nRows = 1; nCols = 1;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
isub = 1;


if 1 % B angle
  hca = h(isub); isub = isub + 1;
  plot3(hca,xTHOR(:),yTHOR(:),zTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),[],AngleBpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = '\Theta_{B}';
  colormap(hca,mirrorcmap);
  hca.CLim = [0 180];
  hca.Title.String = 'Magnetic field normal angle';

  hLegend = [];
  txtLegend = {};
  if plotB
    hold(hca,'on')
    quivB = quiver3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),tsB.z.data(isCrossing,1),'k');
    hold(hca,'off')
    hLegend = [hLegend(:) quivB];
    txtLegend = {txtLegend{:} 'Magnetic field'};
  end
  if plotShockNormal
    hold(hca,'on')
    quivN = quiver3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),shockNormal.x.data(isCrossing,1),shockNormal.y.data(isCrossing,1),shockNormal.z.data(isCrossing,1),'color',[1 0 1]);
    hold(hca,'off')
    hLegend = [hLegend(:) quivN];
    txtLegend = {txtLegend{:} 'Shock normal'};
  end
  if ~isempty(hLegend)
    legend(hLegend,txtLegend,'location','bestoutside')
  end
  hca.FontSize = fontsize;
end

for ii = 1:nRows*nCols
  axis(h(ii),'equal')
  h(ii).XGrid = 'on';
  h(ii).YGrid = 'on';
  h(ii).ZGrid = 'on';
  h(ii).XTick = xtick;
  h(ii).YTick = ytick;
  h(ii).ZTick = xtick;
  h(ii).XLim = xlim;
  h(ii).YLim = ylim;
  h(ii).XLabel.String = 'X';
  h(ii).YLabel.String = 'Y';
  h(ii).FontSize = fontsize;
end
if doPrint, eval(printString); end

%% Figure: 3D: QV QBperp QBpar, for checking the 3D-case
figure_name = 'Q_3D';
figure_position = scrsz; figure_position(3) = figure_position(3)*0.5;
figure('Position',figure_position)

nRows = 2; nCols = 2;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
isub = 1;

if 1 % QV
  hca = h(isub); isub = isub + 1;
  plot3(hca,xTHOR(:),yTHOR(:),zTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),[],QV.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{V}';
  if plotB
    hold(hca,'on')
    quiver3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),tsB.z.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % B angle
  hca = h(isub); isub = isub + 1;
  plot3(hca,xTHOR(:),yTHOR(:),zTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),[],AngleBpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = '\Theta_{B}';
  colormap(hca,mirrorcmap);
  hca.CLim = [0 180];
  hca.Title.String = 'Shock normal angle';

  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  if plotShockNormal
    hold(hca,'on')
    quiver3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),tsB.z.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % QB perp
  hca = h(isub); isub = isub + 1;
  plot3(hca,xTHOR(:),yTHOR(:),zTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),[],QBperp.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{B,\perp}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  if plotShockNormal
    hold(hca,'on')
    quiver3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),tsB.z.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % QB par
  hca = h(isub); isub = isub + 1;
  plot3(hca,xTHOR(:),yTHOR(:),zTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),[],QBpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = sprintf('Q_{B,||}','a');
  if plotB
    hold(hca,'on')
    quiver3(hca,xTHOR(isCrossing),yTHOR(isCrossing),zTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),tsB.z.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end


for ii = 1:nRows*nCols
  axis(h(ii),'equal')
  h(ii).XGrid = 'on';
  h(ii).YGrid = 'on';
  h(ii).XTick = xtick;
  h(ii).YTick = ytick;
  h(ii).XLim = xlim;
  h(ii).YLim = ylim;
  h(ii).XLabel.String = 'X';
  h(ii).YLabel.String = 'Y';
  h(ii).FontSize = fontsize;
end
linkaxes(h,'xy')
if doPrint, eval(printString); end

%% Figure: BAngle, QBperp, QBpar, finite Bz
figure_name = 'theta_QBperp_QBpar';
if 0 % 1x4
  figure_position = scrsz; figure_position(4) = figure_position(4)*0.4;
  figure('Position',figure_position)
  nRows = 1; nCols = 4;
else % 2x2
  figure_position = scrsz;
  figure_position(3) = figure_position(3)*0.6;
  figure_position(4) = figure_position(4)*0.75;
  figure('Position',figure_position)
  nRows = 2; nCols = 2;
end
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
hold(hca,'on')
scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QBperp.data(isCrossing,:),'marker',markerStyle)%,'filled')
hold(hca,'off')
hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
colormap(hca,cmap);
hca.CLim = [0 1];
hca.Title.String = 'Q_{B,\perp}';
if plotB
  hold(hca,'on')
  quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
  hold(hca,'off')
end

hca = h(isub); isub = isub + 1;
plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
hold(hca,'on')
scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QBpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
hold(hca,'off')
hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
colormap(hca,cmap);
hca.CLim = [0 1];
hca.Title.String = sprintf('Q_{B,||}','a');
if plotB
  hold(hca,'on')
  quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
  hold(hca,'off')
end

hca = h(isub); isub = isub + 1;
plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
hold(hca,'on')
scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],AngleBparBz.data(isCrossing,:),'marker',markerStyle)%,'filled')
hold(hca,'off')
hcb = colorbar('peer',hca); hcb.YLabel.String = '\Theta_B';
colormap(hca,mirrorcmap);
hca.Title.String = 'Shock normal angle';
if plotB
  hold(hca,'on')
  quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
  hold(hca,'off')
end

if 1 % Cumlative distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQBperp(end:-1:1)),'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQBpar(end:-1:1)),'histc','g');
  bars.FaceAlpha = 0.5;
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQperp(end:-1:1)),'histc','y');
  bars.FaceAlpha = 0.5;
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQpar(end:-1:1)),'histc','b');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{B,\perp}','Q_{B,||}','Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Cumulative distribution of Q_{B}';
  hca.FontSize = fontsize;
  limQ = 0.8;
  axes(hca)
  infotext{1} = sprintf('Q_{B,perp}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(QBperp(isCrossing).data>limQ)),'%.0f'));
  infotext{2} = sprintf('Q_{B,||}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(QBpar(isCrossing).data>limQ)),'%.0f'));
  infotext{3} = sprintf('Q_{perp}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qperp(isCrossing).data>limQ)),'%.0f'));
  infotext{4} = sprintf('Q_{||}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qpar(isCrossing).data>limQ)),'%.0f'));
  ht=text(1.02,hca.YLim(2)*0.95,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
  ht.HorizontalAlignment = 'left';
  axes(hca)
  infotext = {};
  infotext{1} = sprintf('%s',t1utc(1:10));
  infotext{2} = sprintf('%s',t2utc(1:10));
  ht=text(0.45,hca.YLim(2)*0.99,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
  ht.HorizontalAlignment = 'center';
  ht.FontSize = fontsize;
end

if 0 % Distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ,nDistQperp,'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ,nDistQpar,'histc','g');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Distribution of Q=Q_{R}Q_{V}Q_{B}';
  hca.FontSize = fontsize;
  axis(hca,'square')
end

for ih = 4
  h(ih).YGrid = 'on';
  h(ih).XGrid = 'on';
  axis(h(ii),'square')
end

for ii = 1:3
  axis(h(ii),'equal')
  h(ii).XGrid = 'on';
  h(ii).YGrid = 'on';
  h(ii).XTick = xtick;
  h(ii).YTick = ytick;
  h(ii).XLim = xlim;
  h(ii).YLim = ylim;
  h(ii).XLabel.String = 'X';
  h(ii).YLabel.String = 'Y';
  h(ii).FontSize = fontsize;
end
linkaxes(h(1:3),'xy')

if doPrint, eval(printString); end

%% Figure: All Q: QV QR QBperp QBpar Qperp Qpar
figure_name = 'allQ';
figure_position = scrsz; figure_position(3) = figure_position(3)*0.5;
figure('Position',figure_position)

nRows = 3; nCols = 2;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
isub = 1;

if 1 % QV
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QV.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{V}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % QR
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QR.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{R}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % QB perp
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QBperp.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{B,\perp}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  if plotShockNormal
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),shockNormal.x.data(isCrossing,1),shockNormal.y.data(isCrossing,1),'color',0.5*[1 1 1])
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % QB par
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QBpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = sprintf('Q_{B,||}','a');
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % Qperp
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],Qperp.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{\perp}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % Qpar
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],Qpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = sprintf('Q_{||}','a');
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end

for ii = 1:6
  axis(h(ii),'equal')
  h(ii).XGrid = 'on';
  h(ii).YGrid = 'on';
  h(ii).XTick = xtick;
  h(ii).YTick = ytick;
  h(ii).XLim = xlim;
  h(ii).YLim = ylim;
  h(ii).XLabel.String = 'X';
  h(ii).YLabel.String = 'Y';
  h(ii).FontSize = fontsize;
end
linkaxes(h(1:6),'xy')
if doPrint, eval(printString); end

%% Figure: Summarizing figure
figure_name = 'summary1';
figure_position = scrsz;
figure_position = scrsz; figure_position(4) = figure_position(4)*0.6;
figure('Position',figure_position)

nRows = 2; nCols = 4;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
isub = 1;

if 1 % Distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ,nDistQperp,'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ,nDistQpar,'histc','g');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Distribution of Q=Q_{R}Q_{V}Q_{B}';
  hca.FontSize = fontsize;
  axes(hca)
  infotext{1} = sprintf('%s',t1utc(1:10));
  infotext{2} = sprintf('%s',t2utc(1:10));
  ht=text(0.4,hca.YLim(2)*0.95,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
  ht.HorizontalAlignment = 'center';
  ht.FontSize = fontsize;
end
if 1 % Cumlative distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQperp(end:-1:1)),'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQpar(end:-1:1)),'histc','g');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Cumulative distribution of Q=Q_{R}Q_{V}Q_{B}';
  hca.FontSize = fontsize;
  limQ = 0.8;
  axes(hca)
  infotext{1} = sprintf('Q_{perp}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qperp(isCrossing).data>limQ)),'%.0f'));
  infotext{2} = sprintf('Q_{||}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qpar(isCrossing).data>limQ)),'%.0f'));
  ht=text(0.1,hca.YLim(2)*0.95,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
end

if 1 % QV
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QV.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{V}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % Schock normal angle
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],AngleBparBz0.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = '\Theta_B';
  colormap(hca,mirrorcmap);
  hca.Title.String = {'Angle between shock normal','and B_{IMF}'};
  hca.FontSize = fontsize;
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
end
if 1 % QB perp
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QBperp.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{B,\perp}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % QB par
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],QBpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = sprintf('Q_{B,||}','a');
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1 % Qperp
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],Qperp.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{\perp}';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],Qpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = sprintf('Q_{||}','a');
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end

for ih = 1:2
  h(ih).YGrid = 'on';
  h(ih).XGrid = 'on';
end

for ii = 3:8
  axis(h(ii),'equal')
  h(ii).XGrid = 'on';
  h(ii).YGrid = 'on';
  h(ii).XTick = xtick;
  h(ii).YTick = ytick;
  h(ii).XLim = xlim;
  h(ii).YLim = ylim;
  h(ii).XLabel.String = 'X';
  h(ii).YLabel.String = 'Y';
  h(ii).FontSize = fontsize;
end
linkaxes(h(3:8),'xy')
if doPrint, eval(printString); end

%% Figure: Summarizing figure, distributions of Qx2 and Qpar, Qperp
figure_name = 'summary_dist_Qpar_Qperp';
figure_position = scrsz;
figure_position(3) = figure_position(3)*0.6;
figure_position(4) = figure_position(4)*0.75;
figure('Position',figure_position)


nRows = 2; nCols = 2;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end

isub = 1;

if 1 % Distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ,nDistQperp,'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ,nDistQpar,'histc','g');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Distribution of Q=Q_{R}Q_{V}Q_{B}';
  hca.FontSize = fontsize;
  axes(hca)
  infotext{1} = sprintf('%s',t1utc(1:10));
  infotext{2} = sprintf('%s',t2utc(1:10));
  ht=text(0.4,hca.YLim(2)*0.95,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
  ht.HorizontalAlignment = 'center';
  ht.FontSize = fontsize;
end
if 1 % Cumlative distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQperp(end:-1:1)),'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQpar(end:-1:1)),'histc','g');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Cumulative distribution of Q=Q_{R}Q_{V}Q_{B}';
  hca.FontSize = fontsize;
  limQ = 0.8;
  axes(hca)
  infotext{1} = sprintf('Q_{perp}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qperp(isCrossing).data>limQ)),'%.0f'));
  infotext{2} = sprintf('Q_{||}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qpar(isCrossing).data>limQ)),'%.0f'));
  ht=text(0.1,hca.YLim(2)*0.95,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
end

if 1 % Qperp
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],Qperp.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = 'Q_{\perp} (finite B_Z)';
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,xTHOR(:),yTHOR(:),'color',colorOrbit)
  hold(hca,'on')
  scatter(hca,xTHOR(isCrossing),yTHOR(isCrossing),[],Qpar.data(isCrossing,:),'marker',markerStyle)%,'filled')
  hold(hca,'off')
  hcb = colorbar('peer',hca); hcb.YLabel.String = 'Q';
  colormap(hca,cmap);
  hca.CLim = [0 1];
  hca.Title.String = sprintf('Q_{||} (finite B_Z)','a');
  if plotB
    hold(hca,'on')
    quiver(hca,xTHOR(isCrossing),yTHOR(isCrossing),tsB.x.data(isCrossing,1),tsB.y.data(isCrossing,1),'k')
    hold(hca,'off')
  end
  hca.FontSize = fontsize;
end

for ii = 3:4
  axis(h(ii),'equal')
  h(ii).XGrid = 'on';
  h(ii).YGrid = 'on';
  h(ii).XTick = xtick;
  h(ii).YTick = ytick;
  h(ii).XLim = xlim;
  h(ii).YLim = ylim;
  h(ii).XLabel.String = 'X';
  h(ii).YLabel.String = 'Y';
  h(ii).FontSize = fontsize;
end
linkaxes(h(3:4),'xy')
for ih = 1:2
  h(ih).YGrid = 'on';
  h(ih).XGrid = 'on';
end

%% Figure: Distributions of Qx2
figure_name = 'dist_cdf';
figure_position = scrsz;
figure_position(3) = figure_position(3)*0.6;
figure_position(4) = figure_position(4)*0.45;
figure('Position',figure_position)

nRows = 1; nCols = 2;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
isub = 1;

if 1 % Distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ,nDistQperp,'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ,nDistQpar,'histc','g');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Distribution of Q=Q_{R}Q_{V}Q_{B}';
  hca.FontSize = fontsize;
  axes(hca)
  infotext{1} = sprintf('%s',t1utc(1:10));
  infotext{2} = sprintf('%s',t2utc(1:10));
  ht=text(0.4,hca.YLim(2)*0.95,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
  ht.HorizontalAlignment = 'center';
  ht.FontSize = fontsize;
end
if 1 % Cumlative distribution of quality factors Q
  hca = h(isub); isub = isub + 1;
  %bars = bar(hca,centerAngle,nAngles,1);
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQperp(end:-1:1)),'histc','r');
  bars.FaceAlpha = 0.5;
  hold(hca,'on')
  bars = bar(hca,edgesQ(end:-1:1),cumsum(nDistQpar(end:-1:1)),'histc','g');
  bars.FaceAlpha = 0.5;
  %bar(hca,centerAngle,nAnglesBz0,1)
  hold(hca,'off')
  hleg = legend(hca,'Q_{\perp}','Q_{||}');
  hleg.Box = 'off';
  hca.XLim = [0 1.01];
  hca.XTick = [0:0.1:1];
  hca.XLabel.String = 'Q';
  hca.YLabel.String = 'Counts';
  hca.Title.String = 'Cumulative distribution of Q=Q_{R}Q_{V}Q_{B}';
  hca.FontSize = fontsize;
  limQ = 0.8;
  axes(hca)
  infotext{1} = sprintf('Q_{perp}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qperp(isCrossing).data>limQ)),'%.0f'));
  infotext{2} = sprintf('Q_{||}>%s: %s',num2str(limQ,'%.1f'),num2str(numel(find(Qpar(isCrossing).data>limQ)),'%.0f'));
  ht=text(0.2,hca.YLim(2)*0.95,infotext);
  %ht.Intepreter = 'latex';
  ht.VerticalAlignment = 'top';
end

for ih = 1:2
  h(ih).YGrid = 'on';
  h(ih).XGrid = 'on';
end

if doPrint, eval(printString); end
