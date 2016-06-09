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

%% Download OMNI database data
% 3.3 years (duration of the orbit) of representative solar wind conditions 
% should be chosen
loadDataFromFile = 1;
if loadDataFromFile
  load /Users/Cecilia/MATLAB/irfu-matlab/mission/thor/orbit_coverage/omni_data_20010101-20041231_bsnx_Bxyz_MS.mat
else % download data from omni database, change the years as is appropriate
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

  % Can only seem to download one years data at a time
  clear tintUTC
  tsub = 1;
  c_eval('tintUTC{tsub} = ''200?-01-01T00:00:00/200?-12-31T23:59:00''; tsub = tsub+1;',1:4);

  omni_orig = [];
  tic;
  for iy = 1:numel(tintUTC);  
    tint = irf.tint(tintUTC{iy});
    tmp_omni = irf_get_data_omni(tint,'bsnx,Bx,By,Bz,Ms','omni_min');
    omni_orig = [omni_orig; tmp_omni];
    toc
  end

  % Clean up data
  omni = omni_orig;
  omni(isnan(omni(:,2)),:)=[]; % remove all points that dont have R0 data

  R0 = bsnx(:,2); % RE
  kmR0 = bsnx(:,2)*units.RE*1e-3; % km

  tsBSNX = irf.ts_scalar(irf_time(omni(:,1),'epoch>utc'),omni(:,2)); tsBSNX.units = 'RE'; tsBSNX.name = 'Bowshock nose distance, X'; 
  tsB = irf.ts_vec_xyz(irf_time(omni(:,1),'epoch>utc'),omni(:,3:5)); tsB.units = 'nT'; tsB.name = 'Solar wind magnetic field'; 
  tsM = irf.ts_scalar(irf_time(omni(:,1),'epoch>utc'),omni(:,6)); tsM.units = ''; tsM.name = 'Solar wind Mach number'; 
end


%% Adjust timelines of BSNX and THOR
units = irf_units;
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
fy = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645

%% Check if THOR is inside bowshock or not
if 1
  xTHOR = rTHOR.x.data/units.RE*1e3; % km->RE
  yTHOR = rTHOR.y.data/units.RE*1e3;

  yBS = fy(xTHOR,tsBSNX.data); 
  tsyBS = irf.ts_scalar(tsBSNX.time,yBS);

  allInd = 1:rTHOR.length;
  isInside = find(xTHOR<tsBSNX.data & abs(yBS)>abs(yTHOR));
  isOutside = tocolumn(setdiff(allInd,isInside));

  isInboundCrossing = isInside(find(diff(isInside)>1));
  isOutboundCrossing = isOutside(find(diff(isOutside)>1));
  isCrossing = sort([isOutboundCrossing; isInboundCrossing]);
  nCrossings = numel(isCrossing);
  %nCrossings = numel(isOutboundCrossing)+numel(isInboundCrossing);

  nCrossingsPerYear = nCrossings/((rTHOR.time.stop-rTHOR.time.start)/60/60/24/365);
end

%% Check quality of the bowshock crossing
Rout = 15;
QR = thor_QR(tsBSNX,Rout);
QV = thor_QV(tsBSNX,rTHOR);

% Using Bz=0 is mostly for bugcheck
tsB_Bz0 = tsB; tsB_Bz0.data(:,3) = 0;

[QBparBz0,AngleBparBz0,NormalDirection] = thor_QB(tsBSNX,rTHOR,tsB_Bz0);
[QBpar,AngleBpar] = thor_QB(tsBSNX,rTHOR,tsB);
%[QBpar,AngleBpar,NormalDirection] = thor_QB(tsBSNX(isCrossing),rTHOR(isCrossing),tsB(isCrossing));
[QBperpBz0,AngleBperpBz0] = thor_QB(tsBSNX,rTHOR,tsB_Bz0,'perp');
[QBperp,AngleBperp] = thor_QB(tsBSNX,rTHOR,tsB,'perp');

%AngleBperpBz0.data(,:)


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
