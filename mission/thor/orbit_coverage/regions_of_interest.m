%% THOR Region of Interest Definition
% Defines the regions of interest for THOR and displays the boundaries.
%
%
% Dimensions:
%   
%   Distances in Earth Radii
%   

% If nothing else is specified, use the values from the memo.
if ~exist('useMemoValues','var'); useMemoValues = 1; end

%% CONSTANTS AND DEFINITIONS
% Defines the boundary surfaces for the regions of interest. Boundaries in
% elevation and azimuth faces of the volume enclosed by these surfaces are 
% given by the max,min angles provided and not derived/depicted.

% Sector Angles
azim = pi/180*(0:20:360); % Azimuthal angle [0 360]
% All regions except Shock
elev90 = 0.5*pi-pi/180*(0:5:45); % Angle with respect to the ZY-plame, pi/2 is towards Sun
% Shock region
elev50 = 0.5*pi-pi/180*(0:5:25); % Angle with respect to the ZY-plame, pi/2 is towards Sun

if useMemoValues
    % Magnetosheath Boundaries 
    rInnerMS = [12, 12.01, 12.05, 12.11, 12.20, 12.31, 12.45, 12.63, 12.83, 13.06]; %Define as per Memo
    rOuterMS = [13, 13.02, 13.08, 13.19, 13.35, 13.55, 13.80, 14.11, 14.48, 14.92]; %Define

    % Shock Boundaries
    rInnerShock = [13, 13.02, 13.08, 13.19, 13.35, 13.55];
    rOuterShock = [15, 15.02, 15.08, 15.18, 15.32, 15.51];

    % Foreshock Boundaries
    rInnerFS = [20, 20.01, 20.04, 20.09, 20.17, 20.28, 20.44, 20.65, 20.92, 21.26];
    rOuterFS = [26, 26.02, 26.06, 26.15, 26.19, 26.28, 26.74, 27.04, 27.51, 28.04];

    %  Pristine Solar Wind Boundary (not limited to the outside)
    rInnerPSW = 30;
else % use values from region_coverage.m    
    th_0_45 = round(numel(Th)/2):5:find(Th==45);
    th_0_25 = round(numel(Th)/2):5:find(Th==25);
    rInnerMS = R_innerMS(th_0_45);
    rOuterMS = R_outerMS(th_0_45);
    rInnerShock = R_innerBS(th_0_25);
    rOuterShock = R_outerBS(th_0_25);
    rInnerFS = R_innerFS(th_0_45);
    rOuterFS = R_outerFS(th_0_45);
    rInnerPSW = R_innerSW(round(numel(Th)/2));
end

%% CONVERSIONS
[AZIM90,ELEV90] = meshgrid(azim, elev90);
[AZIM50,ELEV50] = meshgrid(azim, elev50);

% Convert all from spherical to cartesian coordinates
% The output from sph2cart is rotated [x,y,z] -> [y,z,x] so that the 
% azimuthal angle is now in YGSE,ZGSE-plane

% Magnetosheath
[YInnerMS, ZInnerMS, XInnerMS] = sph2cart(AZIM90, ELEV90, repmat(rInnerMS, numel(azim),1)');
[YOuterMS, ZOuterMS, XOuterMS] = sph2cart(AZIM90, ELEV90, repmat(rOuterMS, numel(azim),1)');

% Shock 
[YInnerShock, ZInnerShock, XInnerShock] = sph2cart(AZIM50, ELEV50, repmat(rInnerShock, numel(azim),1)');
[YOuterShock, ZOuterShock, XOuterShock] = sph2cart(AZIM50, ELEV50, repmat(rOuterShock, numel(azim),1)');

% Foreshock
[YInnerFS, ZInnerFS, XInnerFS] = sph2cart(AZIM90, ELEV90, repmat(rInnerFS, numel(azim),1)');
[YOuterFS, ZOuterFS, XOuterFS] = sph2cart(AZIM90, ELEV90, repmat(rOuterFS, numel(azim),1)');

% Pristine Solar Wind
[YInnerPSW, ZInnerPSW, XInnerPSW] = sph2cart(AZIM90, ELEV90, rInnerPSW);

%% PLOT

fullRoI = figure;
fullRoI.Name = 'THOR RoI';

title('THOR Regions Of Interest')
xlabel('X_{GSE}')
ylabel('Y_{GSE}')
zlabel('Z_{GSE}')
hold on

sphere(20);

innerMS = surf(XInnerMS, YInnerMS, ZInnerMS);
    innerMS.FaceColor = [1, .84, 0];
    innerMS.DisplayName = 'Magentosheat Inner Boundary';
    innerMS.FaceAlpha = 0.5;
 
outerMS = surf(XOuterMS, YOuterMS, ZOuterMS);
    outerMS.FaceColor = [1, .84, 0];
    outerMS.DisplayName = 'Magentosheat Outer Boundary';
    outerMS.FaceAlpha = 0.5;
    
innerShock = surf(XInnerShock, YInnerShock, ZInnerShock);
    innerShock.FaceColor = [.69, .08, .18];
    innerShock.DisplayName = 'Shock Inner Boundary';
    innerShock.FaceAlpha = 0.5;
 
outerShock = surf(XOuterShock, YOuterShock, ZOuterShock);
    outerShock.FaceColor = [.69, .08, .18];
    outerShock.DisplayName = 'Shock Outer Boundary';
    outerShock.FaceAlpha = 0.5;
        
innerFS = surf(XInnerFS, YInnerFS, ZInnerFS);
    innerFS.FaceColor = [0, .5, 0];
    innerFS.DisplayName = 'Foreshock Inner Boundary';
    innerFS.FaceAlpha = 0.5;
 
outerFS = surf(XOuterFS, YOuterFS, ZOuterFS);
    outerFS.FaceColor = [0, .5, 0];
    outerFS.DisplayName = 'Foreshock Outer Boundary';
    outerFS.FaceAlpha = 0.5;
       
innerPSW = surf(XInnerPSW, YInnerPSW, ZInnerPSW);
    innerPSW.FaceColor = [0, .45, .74];
    innerPSW.DisplayName = 'Pristine Solar Wind Boundary';
    innerPSW.FaceAlpha = 0.5;
   
plot3([0, rInnerPSW*sind(45)],[0, rInnerPSW*cosd(0)*sind(45)], [0, rInnerPSW*sind(0)*sind(45)] , 'Color','black');
plot3([0, rInnerPSW*sind(45)], [0, rInnerPSW*cosd(90)*sind(45)], [0, rInnerPSW*sind(90)*sind(45)], 'Color','black');
plot3([0, rInnerPSW*sind(45)], [0, rInnerPSW*cosd(180)*sind(45)], [0, rInnerPSW*sind(180)*sind(45)], 'Color','black');
plot3([0, rInnerPSW*sind(45)], [0, rInnerPSW*cosd(270)*sind(45)], [0, rInnerPSW*sind(270)*sind(45)], 'Color','black');

leg = legend('show');
leg.Location = 'eastoutside';

hold off
axis square
axis equal