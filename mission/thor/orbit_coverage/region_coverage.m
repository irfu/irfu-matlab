% Gives the time spent in regions of interest. 
%% get_orbit.m 
% Get orbit coordinates.
get_orbit
orb.t = t;
orb.dt = t(2)-t(1);
orb.r = sqrt(x.^2+y.^2); % in RE
orb.theta = atan2d(y,x); % degrees
% Plot
% hp = polar(orb.theta*pi/180,orb.r); hp.Color = [0.2 0.2 0.2];    

%% Solar wind dynamic pressure and 'alpha'
% Get solar wind dynamic pressure and paramter alpha to calculate 
% magnetopause shape. We find the values that correspond to the chosen
% magnetopause standoff distance.

% Load from omni data.
tint = [irf_time([2013 01 01 01 01 00]) irf_time([2013 12 31 23 59 00])];
disp(['Loading omni data: ' irf_time(tint(1),'epoch>utc_yyyy-mm-ddTHH:MM:SS') ' to ' irf_time(tint(2),'epoch>utc_yyyy-mm-ddTHH:MM:SS')])    
omni_data = irf_get_data_omni(tint,'n,V,Bz','omni_hour');        

%% Look at data
t = omni_data(:,1);
n = omni_data(:,2); % cc
v = omni_data(:,3); v(v>1500) = NaN;% km/s 
Bz = omni_data(:,4); % nT

RE = 6371200; % m
mu0 = 1.2566e-6; % 4*pi*1e-7
B0 = 7.94e22; % Am^2 Earth's dipole moment, B = B0/r^3
mp = 1.6726e-27; % proton mass
m = 0.98*mp+0.02*2*mp; % 98 percent hydrogen, 2 percent helium

Dp = n.*1e6*m.*(v*1e3).^2*1e9; % Dynamic pressure, nPa

% Daily averages
nDays = ceil(diff(tint)/60/60/24);      
av_t = nanmean(reshape([t; nan(24-mod(size(t,1),24),1)],24,nDays),1)';
av_n = nanmean(reshape([n; nan(24-mod(size(n,1),24),1)],24,nDays),1)';
av_v = nanmean(reshape([v; nan(24-mod(size(v,1),24),1)],24,nDays),1)';
av_Bz = nanmean(reshape([Bz; nan(24-mod(size(Bz,1),24),1)],24,nDays),1)';
av_Dp = nanmean(reshape([Dp; nan(24-mod(size(Dp,1),24),1)],24,nDays),1)';

% Magnetopause
% Reference: Shue et al 1998
% Eq.(1) r=rzero*(2/(1+cos(theta)))^alpha
% Eq.(9) rzero=(10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6)
% Eq.(10) alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp))
% Default values: Dp=2nPa, Bz=0nT
% av_Bz = 0; av_Dp = 1.4; % default
rzero=(10.22+1.29*tanh(0.184*(av_Bz+8.14))).*av_Dp.^(-1/6.6); % magnetopause standoff distance
alpha=(0.58-0.007*av_Bz).*(1+0.024*log(av_Dp));

if 0 % plot data
    %% 
    h = irf_plot({[av_t,av_Bz],[av_t,av_n],[av_t,av_v],[av_t,av_Dp],[av_t,rzero],[av_t,alpha]});
    h(1).Title.String = 'Solar wind parameters, daily averages';
    h(1).YLabel.String = 'B_z [nT]';
    h(2).YLabel.String = 'n [cc]';        
    h(3).YLabel.String = 'v [km/s]';
    h(4).YLabel.String = 'D_p [nPa]';
    h(5).YLabel.String = 'R_0 [R_E]';
    h(6).YLabel.String = 'alpha';
    %%
    h = scatter3(av_Dp,av_Bz,rzero,rzero,rzero);
    hc = colorbar;
    h.Parent.XLabel.String = 'D_p [nPa]';
    h.Parent.YLabel.String = 'B_z [nT]';
    h.Parent.ZLabel.String = 'R [R_E]';
end

%% Define all R(theta) that is needed

% Magnetopause: magnetosheath inner boundary
% Reference: Shue et al 1998
% Eq.(1) r=rzero*(2/(1+cos(theta)))^alpha
% Eq.(9) rzero=(10.22+1.29*tanh(0.184*(Bz+8.14)))*Dp^(-1/6.6)
% Eq.(10) alpha=(0.58-0.007*Bz)*(1+0.024*log(Dp))
% Default values: Dp=2nPa, Bz=0nT
% av_Bz = 0; av_Dp = 1.4; % default
% Use Bz = 0; Dp = mean(av_Dp); to get R0, and then scale it up to 12 RE.
Bz = 0; Dp = mean(av_Dp);
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
    hp = polar(orb.theta*pi/180,orb.r);
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
dt = diff(t(1:2)); % how much time one orbit-tick is
TT = NT*orb.dt;
ndays = sum(sum(TT))/60/60/24;
disp(['T = ' num2str(ndays*24) ' hours = ' num2str(ndays) ' days = ' num2str(ndays/365) ' years'])
eval(['time_spent.' regions{iRegion}.region '_days = ndays;'])
regions{iRegion}.DaysSpent = ndays;
end

time_spent

%% Plot the boxes
% Plot orbit
if 1
    [plotx,ploty] = pol2cart(orb.theta*pi/180,orb.r);
    hp = plot(plotx,ploty);
    hca = gca;
    hca.XLabel.String = 'x [Re]';
    hca.YLabel.String = 'y [Re]';
    hca.FontSize = 15;
else
    hp = polar(orb.theta*pi/180,orb.r); 
    hca = gca; 
end
hold(hca,'on');
hp.Color = [0.2 0.2 0.2];
 
% Plot boxes
for iRegion=1:nRegion
    boxTh = [regions{iRegion}.Th regions{iRegion}.Th(end:-1:1) regions{iRegion}.Th(1)];
    boxR =  [regions{iRegion}.R1 regions{iRegion}.R2(end:-1:1) regions{iRegion}.R1(1)];
    hp = polar(boxTh*pi/180,boxR);
    hp.Color = regions{iRegion}.Color;
    hp.LineWidth = 2;
    %legs{iRegion} = regions{iRegion}.region;
    legs{iRegion} = [regions{iRegion}.region ' (' num2str(regions{iRegion}.DaysSpent,'%.1f') ' days)'];
end

axis equal
hca.XLim = [0 80];[0 hca.XLim(2)];
hca.YLim = 30*[-1 1];
legend({'Orbit',legs{:}},'location','northeast')
hold(hca,'off');
if 1
view([180 -90])
legend({'Orbit',legs{:}},'location','northwest')
%hca.YAxisLocation = 'right';
end

