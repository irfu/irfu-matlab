% Calculates the number of bowshock crossing based on the actual position
% of THOR and the bowshock location at the time.

%% Load THOR orbit
%datastore('spice','dir','/Users/Cecilia/calc/SPICE');
orbitKernels = 1;
resampleKernels = 0;
if orbitKernels
  units = irf_units;
  rTHOR = thor_orbit('alt1a.bsp',2*3600);
  if resampleKernels
    newTime = rTHOR.time.start:120:rTHOR.time.stop; % 2 min intervals
    tmpR = rTHOR.resample(newTime);
    rTHOR = tmpR;
  end
else
  get_orbit
  rTHOR = irf.ts_vec_xyz(irf_time(t,'epoch>epochtt'),[x y x*0])
end

%% Download OMNI database data
% 3.3 years (duration of the orbit) of representative solar wind conditions
% should be chosen

if 0
  tStartUTC = '2001-01-01T00:00:00';
  tStart = irf_time(tStartUTC,'utc>epochtt');
  TTHOR = rTHOR.time.stop-rTHOR.time.start;
  TYear = 60*60*24*365;
  %tint = tStart + [0 rTHOR.time.stop-rTHOR.time.start];
  tint = tStart + [0 TTHOR]*0.3;

  % Bowshock nose distance, R0
  omni_bsnx_orig = irf_get_data_omni(tint,'bsnx','omni_min');
  %omni_bsnx_orig = irf_get_data_omni(tint,'bsnx','omni_min');
end

% Can only download one years data at a time
clear tintUTC
tsub = 1;
c_eval('tintUTC{tsub} = ''200?-01-01T00:00:00/200?-12-31T23:59:00''; tsub = tsub+1;',1:4);

bsnx_orig = [];
tic;
for iy = 1:numel(tintUTC)
  tint = irf.tint(tintUTC{iy});
  tmp_bsnx = irf_get_data_omni(tint,'bsnx','omni_min');
  bsnx_orig = [bsnx_orig; tmp_bsnx];
  toc
end

% Clean up data
bsnx = bsnx_orig;
bsnx(isnan(bsnx(:,2)),:)=[];

R0 = bsnx(:,2); % RE
kmR0 = bsnx(:,2)*units.RE*1e-3; % km

tintR0UTC = irf_time(bsnx([1 end],1),'epoch>utc');

%% Adjust timelines of BSNX and THOR
tR0 = irf_time(bsnx(:,1),'epoch>epochTT');
tShift=rTHOR.time.start-tR0.start;
newtR0 = tR0+tShift; % shift the time of bsnx to THOR's time
%newtR0 = newtR0.tlim(rTHOR.time([1 end]));
xBSN = irf.ts_scalar(newtR0,kmR0); xBSN.units = 'km'; xBSN.name = 'Bowshock nose distance';
xBSN = xBSN.tlim(rTHOR.time([1 end]));
rTHOR = rTHOR.resample(xBSN); % upsample orbit times to bsnx's timeline, 1 min

%% Bowshock model
fy = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645

%% Check if THOR is inside bowshock or not
xTHOR = rTHOR.x.data/units.RE*1e3; % km->RE
yTHOR = rTHOR.y.data/units.RE*1e3;

yBS = fy(xTHOR,xBSN.data/units.RE*1e3);
tsyBS = irf.ts_scalar(xBSN.time,yBS);

allInd = 1:rTHOR.length;
isInside = find(xTHOR<xBSN.data/units.RE*1e3 & abs(yBS)>abs(yTHOR));
isOutside = tocolumn(setdiff(allInd,isInside));

isInboundCrossing = isInside(find(diff(isInside)>1));
isOutboundCrossing = isOutside(find(diff(isOutside)>1));
isCrossing = sort([isOutboundCrossing; isInboundCrossing]);
nCrossings = numel(isCrossing);
%nCrossings = numel(isOutboundCrossing)+numel(isInboundCrossing);

nCrossingsPerYear = nCrossings/((rTHOR.time.stop-rTHOR.time.start)/60/60/24/365);

%% Plot results
hca=subplot(1,1,1);
%plot(hca,xTHOR(isOutside),yTHOR(isOutside),'.',...
plot(hca,xTHOR(:),yTHOR(:),'-',...
  xTHOR(isInside),yTHOR(isInside),'.',...
  xTHOR(isCrossing),yTHOR(isCrossing),'.')
axis(hca,'equal'); hca.XLabel.String = 'x (R_E)'; hca.YLabel.String = 'y (R_E)';
tStart = irf.tint(tintUTC{iy});
tOMNI = tStart(1)+[0 xBSN.time.stop-xBSN.time.start];

startUTC = tOMNI.start.utc;
stopUTC = tOMNI.stop.utc;
titleString = {['Bowshock crossings: ' num2str(nCrossings) ', crossings per year: ' num2str(nCrossingsPerYear,'%.0f')],...
  ['(Based on bowshock distance: ' startUTC(1:10) ' to ' stopUTC(1:10) ')']};
hca.Title.String = titleString;
legend('Orbit','Inside bowshock','Bowshock crossings')

hca.FontSize = 18;
hca.YLim = 45*[-1 1];
hca.XLim = 45*[-1 1];

%% Plot results
hca=subplot(1,1,1);
xInside = xTHOR; yInside = yTHOR;
xInside(isOutside) = NaN; yInside(isOutside) = NaN;
xOutside = xTHOR; yOutside = yTHOR;
xOutside(isInside) = NaN; yOutside(isInside) = NaN;

plot(hca,xOutside,yOutside,'-',...
  xInside,yInside,'-',...
  xTHOR(isCrossing),yTHOR(isCrossing),'k.')
axis(hca,'equal'); hca.XLabel.String = 'x (R_E)'; hca.YLabel.String = 'y (R_E)';
tStart = irf.tint(tintUTC{iy});
tOMNI = tStart(1)+[0 xBSN.time.stop-xBSN.time.start];

startUTC = tOMNI.start.utc;
stopUTC = tOMNI.stop.utc;
titleString = {['Bowshock crossings: ' num2str(nCrossings) ', crossings per year: ' num2str(nCrossingsPerYear,'%.0f')],...
  ['(Based on bowshock distance: ' startUTC(1:10) ' to ' stopUTC(1:10) ')']};
hca.Title.String = titleString;
legend('Outside bowshock','Inside bowshock','Bowshock crossings')

hca.FontSize = 18;
hca.YLim = 45*[-1 1];
hca.XLim = 45*[-1 1];
%% Plot results: plot inbound and outbound crossings separately
hca=subplot(1,1,1);
%plot(hca,xTHOR(isOutside),yTHOR(isOutside),'.',...
plot(hca,xTHOR(:),yTHOR(:),'-',...
  xTHOR(isInside),yTHOR(isInside),'.',...
  xTHOR(isOutboundCrossing),yTHOR(isOutboundCrossing),'.',...
  xTHOR(isInboundCrossing),yTHOR(isInboundCrossing),'.')
axis(hca,'equal'); hca.XLabel.String = 'x (R_E)'; hca.YLabel.String = 'y (R_E)';
tStart = irf.tint(tintUTC{iy});
tOMNI = tStart(1)+[0 xBSN.time.stop-xBSN.time.start];

startUTC = tOMNI.start.utc;
stopUTC = tOMNI.stop.utc;
titleString = {['Bowshock crossings: ' num2str(nCrossings) ', crossings per year: ' num2str(nCrossingsPerYear,'%.0f')],...
  ['(Based on bowshock distance: ' startUTC(1:10) '-' stopUTC(1:10) ')']};
hca.Title.String = titleString;
legend('Orbit','Inside bowshock','Inbound crossings','Outbound crossings')

hca.FontSize = 14;

%% Plot results: plot both inbound and outbound crossings separately
hca=subplot(1,1,1);
%plot(hca,xTHOR(isOutside),yTHOR(isOutside),'.',...
plot(hca,xTHOR(:),yTHOR(:),'-',...
  xTHOR(isInside),yTHOR(isInside),'.',...
  xTHOR(isCrossing),yTHOR(isCrossing),'.')
axis(hca,'equal'); hca.XLabel.String = 'x (R_E)'; hca.YLabel.String = 'y (R_E)';
tStart = irf.tint(tintUTC{iy});
tOMNI = tStart(1)+[0 xBSN.time.stop-xBSN.time.start];

startUTC = tOMNI.start.utc;
stopUTC = tOMNI.stop.utc;
titleString = {['Bowshock crossings: ' num2str(nCrossings) ', crossings per year: ' num2str(nCrossingsPerYear,'%.0f')],...
  ['(Based on bowshock distance: ' startUTC(1:10) '-' stopUTC(1:10) ')']};
hca.Title.String = titleString;
legend('Orbit','Inside bowshock','Bowshock crossings')

hca.FontSize = 14;

%% Debugging plots
nRows = 2; nCols = 2;
for ip = 1:nRows*nCols; h(ip) = subplot(nRows,nCols,ip); end
isub = 1;

hca = h(isub); isub = isub+1;
irf_plot(hca,xBSN/units.RE*1e3); irf_zoom(hca,'x',xBSN.time([1 end]))

hca = h(isub); isub = isub+1;
irf_plot(hca,{rTHOR.x/units.RE*1e3,rTHOR.y/units.RE*1e3,tsyBS},'comp');
irf_legend(hca,{'x_{THOR}','y_{THOR}','y_{BS}(x_{THOR},R_0)'},[0.01 0.99])
irf_zoom(hca,'x',rTHOR.time([1 end]))

hca = h(isub); isub = isub+1;
plot(hca,xTHOR(:),yTHOR(:),'-',...
  xTHOR(isOutside),yTHOR(isOutside),'.')
axis(hca,'equal'); hca.XLabel.String = 'x'; hca.YLabel.String = 'y';
hca.Title.String = 'Outside bowshock';

hca = h(isub); isub = isub+1;
plot(hca,xTHOR(:),yTHOR(:),'-',...
  xTHOR(isInside),yTHOR(isInside),'.',...
  xTHOR(isCrossing),yTHOR(isCrossing),'.')
axis(hca,'equal'); hca.XLabel.String = 'x'; hca.YLabel.String = 'y';
hca.Title.String = 'Inside bowshock';

