%% Download data
tint=[irf_time([2013 01 01 1 1 0]) irf_time([2013 12 31 23 59 0])];
omni_data2013 = irf_get_data_omni(tint,'n,V,Bz','omni_hour');
bsnx2013 = irf_get_data_omni(tint,'bsnx','omni_min');
tint=[irf_time([2009 01 01 1 1 0]) irf_time([2009 12 31 23 59 0])];
omni_data2009 = irf_get_data_omni(tint,'n,V,Bz','omni_hour');
bsnx2009 = irf_get_data_omni(tint,'bsnx','omni_min');
tint=[irf_time([2001 01 01 1 1 0]) irf_time([2001 12 31 23 59 0])];
omni_data2001 = irf_get_data_omni(tint,'n,V,Bz','omni_hour');
bsnx2001 = irf_get_data_omni(tint,'bsnx','omni_min');

%% Calculate standoff distances
mp = 1.6726e-27; % proton mass
m = 0.98*mp+0.02*2*mp; % assume solar wind composition is 98 percent hydrogen, 2 percent helium 

bsnx2001(isnan(bsnx2001(:,2)),:) = [];
bsnx2009(isnan(bsnx2009(:,2)),:) = [];
bsnx2013(isnan(bsnx2013(:,2)),:) = [];

% 2001, solar maximum
omni_data = omni_data2001;
t = omni_data(:,1);
n = omni_data(:,2); % cc
v = omni_data(:,3); % km/s
Bz = omni_data(:,4); % nT
Dp = n.*1e6*m.*(v*1e3).^2*1e9; % nPa
rzero2001=(10.22+1.29*tanh(0.184*(Bz+8.14))).*Dp.^(-1/6.6);

% 2009, solar minimum
omni_data = omni_data2009;
t = omni_data(:,1);
n = omni_data(:,2); % cc
v = omni_data(:,3); % km/s
Bz = omni_data(:,4); % nT
Dp = n.*1e6*m.*(v*1e3).^2*1e9; % nPa
rzero2009=(10.22+1.29*tanh(0.184*(Bz+8.14))).*Dp.^(-1/6.6);

% 2013, smaller solar maximum, 2*11 years from 2025
omni_data = omni_data2013;
t = omni_data(:,1);
n = omni_data(:,2); % cc
v = omni_data(:,3); % km/s
Bz = omni_data(:,4); % nT
Dp = n.*1e6*m.*(v*1e3).^2*1e9; % nPa
rzero2013=(10.22+1.29*tanh(0.184*(Bz+8.14))).*Dp.^(-1/6.6);

rzero2001(isnan(rzero2001(:,1)),:) = [];
rzero2009(isnan(rzero2009(:,1)),:) = [];
rzero2013(isnan(rzero2013(:,1)),:) = [];

%% Plot
plotLims = 1; % plot reference distances

hca = subplot(1,1,1); hold(hca,'on');
hpBS = plot(hca,linspace(0,100,size(bsnx2001,1)),sort(bsnx2001(:,2)),...
                linspace(0,100,size(bsnx2009,1)),sort(bsnx2009(:,2)),...
                linspace(0,100,size(bsnx2013,1)),sort(bsnx2013(:,2)));
hpMP = plot(hca,linspace(0,100,size(rzero2001,1)),sort(rzero2001),...
                linspace(0,100,size(rzero2009,1)),sort(rzero2009),...
                linspace(0,100,size(rzero2013,1)),sort(rzero2013));

if plotLims % plot reference distances
    hpLim = plot(hca,[0 100],15*[1 1],'k',...
                     [0 100],13*[1 1],'k',...
                     [0 100],12*[1 1],'k');
    text(4,14,'Bowshock region','verticalalignment','middle','fontsize',14)         
    text(43,12.5,'Magnetosheath region','verticalalignment','middle','fontsize',14)
end
for kk = 1:3 % change linestyle and color to match with BS-color for the same year
    hpMP(kk).Color = hpBS(kk).Color;
    hpMP(kk).LineStyle = '--';
end
hold(hca,'off');

hl = legend(hca,'Bowshock 2001 (OMNI BSNX)','Bowshock 2009 (OMNI BSNX)','Bowshock 2013 (OMNI BSNX)','Magnetopause 2001 ','Magnetopause 2009','Magnetopause 2013','location','northwest');

xlabel('Percentile')
ylabel({'Standoff distance [RE]'}) 
grid on
hca.Box = 'on';
set(gca,'ylim',[7 19])
hca.FontSize = 15;