%% Load THOR orbit
%datastore('spice','dir','/Users/Cecilia/calc/SPICE');
orbitKernels = 1;
resampleKernels = 1;
if orbitKernels
  units = irf_units;
  rTHOR = thor_orbit('new2a.bsp',3600);
  if resampleKernels
    newTime = rTHOR.time.start:120*2:rTHOR.time.stop; % 2 min intervals
    tmpR = rTHOR.resample(newTime);
    rTHOR = tmpR;
  end
else
  get_orbit  %#ok<UNRCH>
  rTHOR = irf.ts_vec_xyz(irf_time(t,'epoch>epochtt'),[x y x*0])
end


%% Different KSRs
% Inside magnetopause, can use theta
Bz = 0; Dp = 2; mpR0 = 12;
R0 = mpR0;
MSP = rTHOR;
thetaTHOR = acos(rTHOR.x.data./sqrt(rTHOR.x.data.^2+rTHOR.y.data.^2+rTHOR.z.data.^2));
radiusTHOR = sqrt(rTHOR.x.data.^2+rTHOR.y.data.^2+rTHOR.z.data.^2)*1e3/units.RE;
alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp));
RR = R0*(2./(1+cos(thetaTHOR))).^alpha;
indInsideMP = find(RR>radiusTHOR); % <------

% Bowshock model
% y is defined from x (GSE)
fy = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
bsR0 = 14; % Bowshock nose distance
R0 = bsR0;
B0 = [-1 1 0]; % Magnetic field, parker spiral
BAng = atan2d(B0(2),B0(1));
theta = linspace(0,180,30); % 'polar' angle, from x_GSE
azimuthal = linspace(0,360,50);
x = R0*cosd(theta);
y = fy(x,R0); % original F/G model adds rstandoff^2=645
r = sqrt(x.^2+y.^2);
theta = atan2d(y,x)*pi/180;

% Transform to 3D
azim = pi/180*azimuthal; % Azimuthal angle [0 360]
elev = pi/2-theta; % Angle with respect to the ZY-plame, pi/2 is towards Sun
[AZ,ELEV] = meshgrid(azim,elev);

% Bowshock surface
[Y,Z,X] = sph2cart(AZ,ELEV,repmat(r,numel(azim),1)');

% Forshock normal
[nX, nY, nZ] = surfnorm(X,Y,Z);
% Find point where magnetic field is tangent to surface
tangentNormal = cross(B0/norm(B0),[0 0 1]); % 45 deg to parker spiral, in xy-plane
tangNormalAngle = acosd(-(nX*tangentNormal(1) + nY*tangentNormal(2) + nZ*tangentNormal(3)));
indTangent = find(tangNormalAngle==min(min(tangNormalAngle)),1,'first');
tX = X(indTangent); % RE
tY = Y(indTangent);
tZ = Z(indTangent);
tAng = tangNormalAngle(indTangent)*180/pi;

% Pick out the different regions
allInd = 1:rTHOR.length;
% Inside bowshock, can use fy(x)
BS = rTHOR;
indNotBS = find(imag(fy(BS.data(:,1)/(units.RE*1e-3),R0))); % Remove all points in front of bowshock nose
indNotBS2 = find(abs(fy(BS.data(:,1)/(units.RE*1e-3),R0))<abs(BS.data(:,2)/(units.RE*1e-3)));
BS.data(indNotBS) = NaN;
BS.data(indNotBS2) = NaN;
indInsideBS = setdiff(allInd,[indNotBS; indNotBS2]);


% Foreshock, outside bowshock but inside solar wind B tangent
% Translate orbit coordinates so that origin is at the tangent point, and
% see what points end up within the foreschock
rTHOR_ = rTHOR-[tX tY tZ]*units.RE*1e-3;
% 0 deg towards Sun (x_GSE)
angleXGSE = atan2d(rTHOR_.data(:,2),rTHOR_.data(:,1));
FS = rTHOR;
indNotFS = find(angleXGSE>BAng-180);
indNotFS2 = find(angleXGSE<-BAng-180);
indInsideFS = setdiff([allInd],[indNotFS;indNotFS2]);
FS.data(indNotFS) = NaN;
FS.data(indNotFS2) = NaN;
FS.data(indInsideBS) = NaN;

indInsideFS = setdiff(allInd,[tocolumn(indInsideBS);tocolumn(indNotFS);tocolumn(indNotFS2)]);

% Collect the various indices and make vectors
% Magnetosphere
indMSP = indInsideMP;
rMSP = rTHOR; rMSP.data(setdiff(allInd,indMSP),:) = NaN;
% Magnetosheath
indMSH = setdiff(indInsideBS,indInsideMP);
rMSH = rTHOR; rMSH.data(setdiff(allInd,indMSH),:) = NaN;
% Foreshock
indFS = setdiff(indInsideFS,indInsideBS);
rFS = rTHOR; rFS.data(setdiff(allInd,indFS),:) = NaN;
% Solar wind, outside bowshock and foreshock
indSW = setdiff(allInd,[tocolumn(indMSP);tocolumn(indMSH);tocolumn(indFS)]);
rSW = rTHOR; rSW.data(setdiff(allInd,indSW)) = NaN;


% Calculate the time spent in each region
dt = rTHOR.time(2)-rTHOR.time(1); % s
tTotal = rTHOR.time.stop-rTHOR.time.start; % s
nMSP = sum(~isnan(rMSP.data(:,1))); tMSP = nMSP*dt;
nMSH = sum(~isnan(rMSH.data(:,1))); tMSH = nMSH*dt;
nFS = sum(~isnan(rFS.data(:,1))); tFS = nFS*dt;
nSW = sum(~isnan(rSW.data(:,1))); tSW = nSW*dt;

% Checking total time spent, need to take away one dt due to endpoints
tDiff = tTotal - ((tSW+tFS+tMSH+tMSP)-dt);


%% Define all R(theta) that is needed

% Magnetopause: magnetosheath inner boundary
% Reference: Shue et al 1998
% Eq.(1) r=rzero*(2/(1+cos(theta)))^alpha
% Eq.(9) rzero=(10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6)
% Eq.(10) alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp))
% Default values: Dp=2nPa, Bz=0nT
% av_Bz = 0; av_Dp = 1.4; % default
% Use Bz = 0; Dp = mean(av_Dp); to get R0, and then scale it up to 12 RE.

Bz = 0; Dp = 2;
rzero = (10.22+1.29*tanh(0.184*(Bz+8.14))).*Dp.^(-1/6.6);
alpha = (0.58-0.007*Bz).*(1+0.024*log(Dp));
r_fun = @(R0,theta,alpha) R0.*(2./(1+cosd(theta))).^alpha;
Th_innerMS = -45:1:45;
R_MP = r_fun(rzero,Th_innerMS,alpha);
R_innerMS = R_MP*12/R_MP(46); % Rescale to R0 = 12 RE

% Magnetosheath outer boundary, bowshock inner boundary
R0 = 13;
theta = 0:45;
x = R0*cosd(theta);
y = sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
r = sqrt(x.^2+y.^2);
theta = atan2d(y,x)*pi/180;
R_innerBS = [r(end:-1:2) r]; R_outerMS = R_innerBS;
Th_innerBS = [-theta(end:-1:2) theta]; Th_outerMS = Th_innerBS;

% Bowshock outer boundary
R0 = 15;
theta = 0:45;
x = R0*cosd(theta);
y = sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
r = sqrt(x.^2+y.^2);
theta = atan2d(y,x)*pi/180;
R_outerBS = [r(end:-1:2) r]; % outer radius
Th_outerBS = [-theta(end:-1:2) theta];

% Forshock inner boundary
R0 = 20;
theta = 0:45;
x = R0*cosd(theta);
y = sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645
r = sqrt(x.^2+y.^2);
theta = atan2d(y,x)*pi/180;
R_innerFS = [r(end:-1:2) r];
Th_innerFS = [-theta(end:-1:2) theta];

% Foreshock outer boundary
R0_outerFS = 26;
R_outerFS = R_innerFS*R0_outerFS/R0; % Rescale from 20 to 26
Th_outerFS = Th_innerFS;

% interpolate values to given thetas
Th = -45:1:45; nTh = numel(Th); cTh = round(nTh/2);
R_outerMS = spline(Th_outerMS'*180/pi,[R_innerBS'],Th);
R_innerBS = spline(Th_innerBS'*180/pi,[R_innerBS'],Th);
R_outerBS = spline(Th_innerBS'*180/pi,[R_outerBS'],Th);
R_innerFS = spline(Th_innerFS'*180/pi,[R_innerFS'],Th);
R_outerFS = spline(Th_outerFS'*180/pi,[R_outerFS'],Th);

% Pristine solar wind inner radius
R_innerSW = ones(1,nTh)*30;
R_outerSW = ones(1,nTh)*70;

if 0 % plot the boundaries
  %%
  hp = polar(orb.theta*pi/180,orb.r); %#ok<UNRCH>
  hp.Color = [0.2 0.2 0.2];
  hold on
  polar(Th*pi/180, R_innerMS);
  polar(Th*pi/180, R_outerMS);
  polar(Th*pi/180, R_innerBS);
  polar(Th*pi/180, R_outerBS);
  polar(Th*pi/180, R_innerFS);
  polar(Th*pi/180, R_outerFS);
  polar(Th*pi/180, R_innerSW);
  polar(Th*pi/180, R_outerSW);
  hold off
  legend('Orbit','Inner Magnetosheath','Outer Magnetosheath','Inner Bowshock','Outer Bowshock','Inner Foreshock','Outer Foreshock','Solar wind','Outer limit of solar wind')
end

% Display values for R
thshow = cTh:5:nTh;
disp(' ')
disp(['   Theta     R_innerMS R_outerMS R_innerBS R_outerBS R_innerFS R_outerFS R_innerSW'])
disp('   -----------------------------------------------------------------------------')
disp([Th(thshow)' R_innerMS(thshow)' R_outerMS(thshow)' R_innerBS(thshow)' R_outerBS(thshow)' R_innerFS(thshow)' R_outerFS(thshow)' R_innerSW(thshow)'])
disp('   Obs! R_innerBS and R_outerBS for Theta>25 are not used.')
%% Define edges of boxes
% Magnetosheath, Th spans 90 deg
iMS = 1:1:91;
Magnetosheath.region = 'Magnetosheath';
Magnetosheath.R1 = R_innerMS(iMS);
Magnetosheath.R2 = R_outerMS(iMS);
Magnetosheath.Th = Th(iMS);
Magnetosheath.Color = [1 0.8 0.0];
% Bowshock, Th spans 40 deg
iBS = 21:1:71;
%iBS = 1:91;
Bowshock.region = 'Bowshock';
Bowshock.R1 = R_innerBS(iBS); % -25:25
Bowshock.R2 = R_outerBS(iBS);
Bowshock.Th = Th(iBS);
Bowshock.Color = [1 0 0.0];
% Foreshock, Th spans 90 deg
Foreshock.region = 'Foreshock';
Foreshock.R1 = R_innerFS;
Foreshock.R2 = R_outerFS;
Foreshock.Th = Th;
Foreshock.Color = [0 0.8 0.0];
% Solar wind, Th spans 90 deg
Solarwind.region = 'Solarwind';
Solarwind.R1 = R_innerSW;
Solarwind.R2 = R_outerSW;
Solarwind.Th = Th;
Solarwind.Color = [0 0 1];
regions = {Magnetosheath,Bowshock,Foreshock,Solarwind};

%% Do binning
orb.t = rTHOR.time.epochUnix;
orb.dt = orb.t(2)-orb.t(1);
x = rTHOR.x.data*1e3/units.RE;
y = rTHOR.y.data*1e3/units.RE;
orb.r = sqrt(x.^2+y.^2); % in RE
orb.theta = atan2d(y,x); % degrees

nRegion = numel(regions);
for iRegion = 1:nRegion
  disp(['------ ' regions{iRegion}.region])
  disp(['Theta = 0: R1 = ' num2str(regions{iRegion}.R1(round(0.5*numel(regions{iRegion}.Th)))) '  R2 = ' num2str(regions{iRegion}.R2(round(0.5*numel(regions{iRegion}.Th))))])
  edgesTh = regions{iRegion}.Th';
  edgesR = [regions{iRegion}.R1' regions{iRegion}.R2'];
  centerTh = edgesTh(1:end-1)+0.5*(edgesTh(2)-edgesTh(1));
  centerR = [spline(edgesTh,edgesR(:,1),centerTh) spline(edgesTh,edgesR(:,2),centerTh)];
  nBinsTh = numel(centerTh);
  NT=0;
  % How much time is spent in each bin
  for kk = 1:nBinsTh
    [nt,edges,mid,loc] = histcn([orb.r(:) orb.theta(:)],centerR(kk,:),edgesTh([kk kk+1],:));
    nts(kk) = nt;
    NT = NT+nt;
  end
  dt = diff(orb.t(1:2)); % how much time one orbit-tick is
  TT = NT*orb.dt;
  ndays = sum(sum(TT))/60/60/24;
  disp(['T = ' num2str(ndays*24) ' hours = ' num2str(ndays) ' days = ' num2str(ndays/365) ' years'])
  eval(['time_spent.' regions{iRegion}.region '_days = ndays;'])
  regions{iRegion}.DaysSpent = ndays;
end

time_spent


%% Plot KSRs
h = subplot(1,1,1);

% Earth
sphere(h,20);
hold(h,'on')

colors = get(h,'colororder');

% Point where the solar wind magnetic field is tangent to the magnetopause
% h_tp = plot3(tX,tY,tZ,'ko');
% h_tp.MarkerSize = 10;
% h_tp.MarkerFaceColor = [0 0 0];

% Plot part of the orbit that is within the magnetopause
h_sw = plot3(h,rSW.data(:,1)/(units.RE*1e-3),rSW.data(:,2)/(units.RE*1e-3),rSW.data(:,3)/(units.RE*1e-3)); h_sw.Color = colors(1,:);
h_fs = plot3(h,rFS.data(:,1)/(units.RE*1e-3),rFS.data(:,2)/(units.RE*1e-3),rFS.data(:,3)/(units.RE*1e-3)); h_fs.Color = colors(5,:);
h_msp = plot3(h,rMSP.data(:,1)/(units.RE*1e-3),rMSP.data(:,2)/(units.RE*1e-3),rMSP.data(:,3)/(units.RE*1e-3)); h_msp.Color = [0.8 0.8 0.8];
h_msh = plot3(h,rMSH.data(:,1)/(units.RE*1e-3),rMSH.data(:,2)/(units.RE*1e-3),rMSH.data(:,3)/(units.RE*1e-3)); h_msh.Color = colors(3,:);

[xMP,yMP] = boundary(12,'mp');
[xBS,yBS] = boundary(14,'bs');
h_bs = plot3(h,xBS,yBS,xBS*0+13,'linewidth',3); h_bs.Color = colors(2,:);
h_mp = plot3(h,xMP,yMP,xMP*0+13,'linewidth',3); h_mp.Color = colors(6,:);

view([0 0 1])
hold(h,'off')

% Add TPRs
for iRegion=1:nRegion
  boxTh = [regions{iRegion}.Th regions{iRegion}.Th(end:-1:1) regions{iRegion}.Th(1)];
  boxR =  [regions{iRegion}.R1 regions{iRegion}.R2(end:-1:1) regions{iRegion}.R1(1)];
  if 0
    [plotx,ploty] = pol2cart(boxTh*pi/180,boxR); %#ok<UNRCH>
    hp = patch(plotx,ploty,'k');
    hp.FaceColor = regions{iRegion}.Color;
    hp.EdgeColor = regions{iRegion}.Color;
    hp.FaceAlpha = 0.2;
  elseif 0 %#ok<IFCDUP>
    hp = polar(boxTh*pi/180,boxR); %#ok<UNRCH>
    hp.Color = regions{iRegion}.Color;
    hp.LineWidth = 2;
  else
    [plotx,ploty] = pol2cart(boxTh*pi/180,boxR);
    hold(h,'on')
    hp_ = plot3(plotx,ploty,plotx*0+13,'k','linewidth',3);
    hp(iRegion) = plot3(plotx,ploty,plotx*0+13,'k-','linewidth',1);
    hp(iRegion).Color = regions{iRegion}.Color;
    hold(h,'off')
  end
  %legs{iRegion} = regions{iRegion}.region;
  legs{iRegion} = [regions{iRegion}.region ' (' num2str(regions{iRegion}.DaysSpent,'%.0f') ' days)'];
end

title(h,{'Key Science Regions (KSR) and Top Priority Regions (TPR)', ['Bowshock nose at ' num2str(bsR0,'%.1f') 'R_E'],['Magnetopause at ' num2str(mpR0,'%.1f') 'R_E']})
title(h,{'Key Science Regions (KSR) and Top Priority Regions (TPR)'})
title(h,{'Key Science Regions (KSRs)', 'and', 'Top Priority Regions (TPRs)'})
xlabel(h,'X_{GSE}')
ylabel(h,'Y_{GSE}')
zlabel(h,'Z_{GSE}')
h.FontSize = 14;

axis equal;
box(h,'on')
h.YLim = [-30 30];
h.XLim = [-10 50];
if 0
  legend([h_sw h_fs h_msh],{['KSR: Solar wind: ' num2str(tSW/60/60/24,'%.0f') ' days'],...
    ['KSR: Foreschock: ' num2str(tFS/60/60/24,'%.0f') ' days'],...
    ['KSR: Magnetosheath: ' num2str(tMSH/60/60/24,'%.0f') ' days']},...
    'location','northeastoutside') %#ok<UNRCH>
end
