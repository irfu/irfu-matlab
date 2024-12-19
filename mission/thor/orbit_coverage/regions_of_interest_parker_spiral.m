% Loads the THOR orbit and calculates the time spent in each of the three
% regions: (1) inside bowshock, (2) outside bowshock but inside foreshock,
% (3) solar wind, outside bowshock and foreshock.

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

%% Bow shock model
% y is defined from x (GSE)
fy = @(x,R0) sqrt(0.04*(x-R0).^2-45.3*(x-R0)); % original F/G model adds rstandoff^2=645

R0 = 14; % Bowshock nose distance
B0 = [-1 1 0]; % Magnetic field, parker spiral, can be changed to some extent
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

indBS = setdiff(allInd,[indNotBS; indNotBS2]);

% Foreshock, outside bowshock but inside solar wind B tangent
% Translate orbit coordinates so that origin is at the tangent point, and
% see what points end up within the foreschock
rTHOR_ = rTHOR-[tX tY tZ]*units.RE*1e-3;
% 0 deg towards Sun (x_GSE)
angleXGSE = atan2d(rTHOR_.data(:,2),rTHOR_.data(:,1));
FS = rTHOR;
indNotFS = find(angleXGSE>BAng-180);
indNotFS2 = find(angleXGSE<-BAng-180);
FS.data(indNotFS) = NaN;
FS.data(indNotFS2) = NaN;
FS.data(indBS) = NaN;

indFS = setdiff(allInd,[tocolumn(indBS);tocolumn(indNotFS);tocolumn(indNotFS2)]);

% Solar wind, outside bowshock and foreshock
SW = rTHOR;
indSW = setdiff(allInd,[tocolumn(indBS);tocolumn(indFS)]);
indNotSW = setdiff(allInd,indSW);
SW.data(indNotSW) = NaN;

% Calculate the time spent in each region
dt = rTHOR.time(2)-rTHOR.time(1); % s
tTotal = rTHOR.time.stop-rTHOR.time.start; % s
nBS = sum(~isnan(BS.data(:,1))); tBS = nBS*dt;
nFS = sum(~isnan(FS.data(:,1))); tFS = nFS*dt;
nSW = sum(~isnan(SW.data(:,1))); tSW = nSW*dt;

% Checking total time spent
% need to take away one dt due to endpoints
tDiff = tTotal - ((tSW+tFS+tBS)-dt);

%% Plot orbit in the different regions
h = subplot(1,1,1);

% Earth
sphere(h,20);
hold(h,'on')

% Bowshock
bs = surf(h,X,Y,Z);
bs.FaceColor = [1 1 1];
bs.DisplayName = 'Magentosheat Inner Boundary';
bs.FaceAlpha = 0.5;

axis(h,'square')
axis(h,'equal')

% Point where the solar wind magnetic field is tangent to the magnetopause
h_tp = plot3(tX,tY,tZ,'ko');
h_tp.MarkerSize = 10;
h_tp.MarkerFaceColor = [0 0 0];

% Plot part of the orbit that is within the magnetopause
h_sw = plot3(h,SW.data(:,1)/(units.RE*1e-3),SW.data(:,2)/(units.RE*1e-3),SW.data(:,3)/(units.RE*1e-3));
h_fs = plot3(h,FS.data(:,1)/(units.RE*1e-3),FS.data(:,2)/(units.RE*1e-3),FS.data(:,3)/(units.RE*1e-3));
h_bs = plot3(h,BS.data(:,1)/(units.RE*1e-3),BS.data(:,2)/(units.RE*1e-3),BS.data(:,3)/(units.RE*1e-3));



view([0 0 1])
hold(h,'off')

title(h,['Time spent in regions of interest: Bowshock nose at ' num2str(R0,'%.1f') 'R_E'])
xlabel(h,'X_{GSE}')
ylabel(h,'Y_{GSE}')
zlabel(h,'Z_{GSE}')

legend([h_sw h_fs h_bs],{['Solar wind: ' num2str(tSW/60/60/24,'%.1f') ' days'],...
  ['Foreschock: ' num2str(tFS/60/60/24,'%.1f') ' days'],...
  ['Inside bowshock: ' num2str(tBS/60/60/24,'%.1f') ' days']},...
  'location','northeastoutside')


%% Plot the different regions individually
h = subplot(2,2,1);
h_sw = plot3(h,SW.data(:,1)/(units.RE*1e-3),SW.data(:,2)/(units.RE*1e-3),SW.data(:,3)/(units.RE*1e-3),...
  FS.data(:,1)/(units.RE*1e-3),FS.data(:,2)/(units.RE*1e-3),FS.data(:,3)/(units.RE*1e-3),...
  BS.data(:,1)/(units.RE*1e-3),BS.data(:,2)/(units.RE*1e-3),BS.data(:,3)/(units.RE*1e-3));
view([0 0 1])
title(h,'Spacecraft trajectory')
h = subplot(2,2,2);
h_sw = plot3(h,SW.data(:,1)/(units.RE*1e-3),SW.data(:,2)/(units.RE*1e-3),SW.data(:,3)/(units.RE*1e-3));
view([0 0 1])
title(h,'Spacecraft trajectory: Solar wind')
h = subplot(2,2,3);
h_fs = plot3(h,FS.data(:,1)/(units.RE*1e-3),FS.data(:,2)/(units.RE*1e-3),FS.data(:,3)/(units.RE*1e-3));
view([0 0 1])
title(h,'Spacecraft trajectory: Foreshock')
h = subplot(2,2,4);
h_bs = plot3(h,BS.data(:,1)/(units.RE*1e-3),BS.data(:,2)/(units.RE*1e-3),BS.data(:,3)/(units.RE*1e-3));
view([0 0 1])
title(h,'Spacecraft trajectory: Within bowshock')


