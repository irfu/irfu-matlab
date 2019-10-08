%% Plot Multiple Histograms
% Generate two vectors of random numbers and plot a histogram for each
%% THOR probability histogram frequencies - SW,MSH 1) get data
%Using Jan's database
dirTHORGoogleDrive = '/Users/andris/Google Drive/Job/Projects/THOR';

cd([dirTHORGoogleDrive '/Figures/matlab/'])

units = irf_units;
% get MSH data
load('../../Plasma parameters/databases/msht_stats_cl1_nov06jun08_cis_ok.mat')
n_msh = pars.dens;
b_msh = pars.bmag;
tperp_msh = pars.tperp*1e6*units.kB/units.e/5; % ion density is given in MK, divide by 5 to get electron temp
tpar_msh = pars.tpar;
T_msh = tperp_msh;
v_vec_msh = pars.vel_gse;
v_msh = sqrt(sum(v_vec_msh.^2,1));
v_te_msh = sqrt(T_msh*units.eV*2/units.me)/1000; % km/s

f_pe_msh = 9000*sqrt(n_msh); % electron plasma frequency
f_ce_msh = 28*b_msh; % electron cyclotron frequency
f_ci_msh = f_ce_msh/1836; % proton cyclotron frequency

de_msh = 7.43*sqrt(T_msh./n_msh)/1000; % electron Debye length, km
cwpi_msh = 230 ./ sqrt(n_msh); % proton inertial length
f_cwpi_msh = v_msh./cwpi_msh; % doppler shifted proton inertial length
f_cwpe_msh = f_cwpi_msh*43; % doppler shifted electron inertial length
f_de_msh = v_te_msh./de_msh; % should be d oppler shifted by electron thermal velocity

roe_msh = 3.4*sqrt(tperp_msh)./b_msh; % electron cyclotron radius
f_roe_msh = v_msh./roe_msh; % doppler shifter electron cyclotron radius
% doppler shifter ion cyclotron radius

% get solar wind data
tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 12 31 23 59 0])];
load ff.mat
load([dirTHORGoogleDrive '/Figures/omni_BsnxBxByBzMsVNT_2006.mat'])

%ff = irf_get_data_omni(tint,'n,v,b');
n_sw = ff(:,2);
v_sw = ff(:,3);
b_sw = ff(:,4);
T_sw = tmp_omni(:,[1 9]); T_sw(isnan(T_sw(:,2)),:) = []; T_sw = irf_resamp(T_sw,ff); T_sw = T_sw(:,2)*units.kB/units.e;
T_sw = T_sw*2; % electron temp is ~2 times larger than ions in sw

v_te_sw = sqrt(T_sw*units.eV*2/units.me)/1000; % km/s

cwpi_sw = 230 ./ sqrt(n_sw);
de_sw = 7.43*sqrt(T_sw./n_sw)/1000; % electron Debye length, km
f_cwpi_sw = v_sw./cwpi_sw;
f_cwpe_sw = f_cwpi_sw*43;
f_pe_sw = 9000*sqrt(n_sw);
f_ce_sw = 28*b_sw;
f_ce_sw = 28*b_sw; % electron cyclotron frequency
f_ci_sw = f_ce_sw/1836; % proton cyclotron frequency
f_de_sw = v_te_sw./de_sw;


% THOR histogram 2 plot
figure(10);
h=irf_plot(2);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);

color_ions  = [1 0.5 0.5];
color_e     = [0.5 0.5 0.9];
color_fp    = [0.7 0.7 0.7];
color_de    = [0.6 0.8 0.5];

xmin = -0.0;
xmax = 5.3;
ymin = 0;
ymax = 0.199;


% Magnetosheath
hca = subplot(2,1,1); h(1) = hca;

h1 = histogram(hca,log10(f_cwpi_msh));
h1.Normalization = 'probability';
h1.BinWidth = 0.05;
h1.FaceColor = color_ions;

hold on;
h2 = histogram(hca,log10(f_cwpe_msh));
h2.Normalization = 'probability';
h2.BinWidth = 0.05;
h2.FaceColor = color_e;

h3 = histogram(hca,log10(f_pe_msh));
h3.Normalization = 'probability';
h3.BinWidth = 0.05;
h3.FaceColor = color_fp;

h4 = histogram(hca,log10(f_ce_msh));
h4.Normalization = 'probability';
h4.BinWidth = 0.05;
h4.FaceColor = color_fp;

h5 = histogram(hca,log10(f_ci_msh));
h5.Normalization = 'probability';
h5.BinWidth = 0.05;
h5.FaceColor = color_fp;

% h6 = histogram(hca,log10(f_de_msh));
% h6.Normalization = 'probability';
% h6.BinWidth = 0.05;
% h6.FaceColor = color_de;

ylabel('probability')

xlim([xmin xmax])
ylim([ymin ymax])

grid on;
set(hca,'XMinorGrid','off')
set(hca,'YMinorGrid','off')

%text(5.1,0.15,{'Debye','scales'},'color',color_de.^3);
fontsize = 14;
text(0.1,0.17,{'H+','scales'},'color',color_ions.^3,'fontsize',fontsize);
text(1.6,0.15,{'e-','scales'},'color',color_e.^3,'fontsize',fontsize);
text(4.0,0.15,{'plasma','frequency'},'color',color_fp.^3,'fontsize',fontsize);
text(2.9,0.07,{'e- cyclotron','frequency'},'color',color_fp.^3,'fontsize',fontsize);
if xmin<0
  text(-0.9,0.14,{'H+ cyclotron','frequency'},'color',color_fp.^3,'fontsize',fontsize);
end

irf_legend('Magnetosheath',[0.85 0.95]);

% Pristine solar wind
hca = subplot(2,1,2); h(2) = hca;

h1 = histogram(hca,log10(f_cwpi_sw));
h1.Normalization = 'probability';
h1.BinWidth = 0.05;
h1.FaceColor = color_ions;

hold on;

h2 = histogram(hca,log10(f_cwpe_sw));
h2.Normalization = 'probability';
h2.BinWidth = 0.05;
h2.FaceColor = color_e;

h3 = histogram(hca,log10(f_pe_sw));
h3.Normalization = 'probability';
h3.BinWidth = 0.05;
h3.FaceColor = color_fp;

h4 = histogram(hca,log10(f_ce_sw));
h4.Normalization = 'probability';
h4.BinWidth = 0.05;
h4.FaceColor = color_fp;

h5 = histogram(hca,log10(f_ci_sw));
h5.Normalization = 'probability';
h5.BinWidth = 0.05;
h5.FaceColor = color_fp;

% h6 = histogram(hca,log10(f_de_sw));
% h6.Normalization = 'probability';
% h6.BinWidth = 0.05;
% h6.FaceColor = color_de;

ylabel('probability')

xlim([xmin xmax])
ylim([ymin ymax])

grid on;
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

irf_legend('Pristine solar wind',[0.85 0.95]);

% text(1.2,0.15,{'Debye','scales'},'color',color_de.^3);
% text(xmin+0.8,0.15,{'H+','scales'},'color',color_ions.^3);
% text(xmin+2.5,0.15,{'e-','scales'},'color',color_e.^3);
% text(xmin+4,0.15,{'plasma','frequency'},'color',color_fp.^3);
% text(xmin+2.7,0.03,{'e- cyclotron','frequency'},'color',color_fp.^3);

xlabel('(V_{flow} / L) or frequency [Hz]')
h(1).XTickLabel = '';
h(2).XTick = (fix(xmin):fix(xmax));
h(1).XTick = h(2).XTick;
h(2).XTickLabel = arrayfun(@(x) {['10^{' num2str(x) '}']}, h(2).XTick);

h(2).Position(2) = h(2).Position(2)+0.11;

% Add instrument sampling
linewidth = 2;
Ts = 0.050;
%plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',1,'color','r');
plot(h(2),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','r');
text(log10(1./Ts),ymax-0.001,'Ts=50ms','rotation',90,...
	'verticalalignment','bottom','horizontalalignment','right',...
	'fontsize',12,'Parent', h(2));

Ts = 0.150;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','r');
%plot(h(2),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',1,'color','r');
text(log10(1./Ts),ymax-0.001,'Ts=150ms','rotation',90,...
	'verticalalignment','bottom','horizontalalignment','right',...
	'fontsize',12,'Parent', h(1));


Ts = 0.005;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','b');
plot(h(2),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','b');
text(log10(1./Ts),ymax-0.001,'Ts=5ms','rotation',90,...
	'verticalalignment','bottom','horizontalalignment','right',...
	'fontsize',12,'Parent', h(1));


Ts = 1e-5;
%plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',1,'color','r');
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','k');
plot(h(2),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','k');
% text(log10(1./Ts),ymax-0.00,['Ts=10' char(181) 's'],'rotation',90,...
% 	'verticalalignment','bottom','horizontalalignment','right',...
%   'fontsize',12,'Parent', h(1));


for Ts = [0.15 0.005 1e-5]
  text(log10(1./Ts),ymax-0.001,'requirement','rotation',90,...
	'verticalalignment','top','horizontalalignment','right',...
	'fontsize',12,'Parent', h(1));
end
for Ts = [0.05]
  text(log10(1./Ts),ymax-0.001,'requirement','rotation',90,...
	'verticalalignment','top','horizontalalignment','right',...
	'fontsize',12,'Parent', h(2));
end
% axis(h(2))
% Ts = 0.030;
% plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',1,'color','b','linestyle',':');
% plot(h(2),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',1,'color','b','linestyle',':');
% text(log10(1./Ts),ymax-0.01,'Ts=30ms','rotation',90,...
% 	'verticalalignment','bottom','horizontalalignment','right',...
% 	'fontsize',12);

%% Single panel with only 3 peaks

% THOR histogram 2 plot
figure('Position',[10 10 800 300]);
set(gcf,'defaultAxesFontSize',18);
set(gcf,'defaultTextFontSize',18);

color_ions  = [1 0.5 0.5];
color_e     = [0.5 0.5 0.9];
color_fp    = [0.7 0.7 0.7];
color_de    = [0.6 0.8 0.5];

xmin = -0.8;
xmax = 5.3;
ymin = 0;
ymax = 0.179;

% Magnetosheath
hca = subplot(1,1,1);
h(1) = hca;

h1 = histogram(hca,log10(f_cwpi_msh));
h1.Normalization = 'probability';
h1.BinWidth = 0.05;
h1.FaceColor = color_ions;

hold on;
h2 = histogram(hca,log10(f_cwpe_msh));
h2.Normalization = 'probability';
h2.BinWidth = 0.05;
h2.FaceColor = color_e;

h3 = histogram(hca,log10(f_pe_msh));
h3.Normalization = 'probability';
h3.BinWidth = 0.05;
h3.FaceColor = color_fp;

ylabel('probability')

xlim([xmin xmax])
ylim([ymin ymax])

grid on;
set(hca,'XMinorGrid','off')
set(hca,'YMinorGrid','off')

%text(5.1,0.15,{'Debye','scales'},'color',color_de.^3);
fontsize = 18;
text(0.05,0.16,{'protons'},'color',color_ions.^3,'fontsize',fontsize);
text(1.52,0.16,{'electrons'},'color',color_e.^3,'fontsize',fontsize);
text(4.1,0.15,{'plasma','frequency'},'color',color_fp.^3,'fontsize',fontsize);

irf_legend('Magnetosheath',[0.6 1.03]);


xlabel('          V_{flow} / L_{kinetic} [1/s]                                            frequency [Hz]')

h(1).XTickLabel = '';
h(1).XTick = (fix(xmin):fix(xmax));
h(1).XTick = h(1).XTick;
h(1).XTickLabel = arrayfun(@(x) {['10^{' num2str(x) '}']}, h(1).XTick);

% Add instrument sampling
linewidth = 2;

Ts = 0.150;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','r');

Ts = 0.005;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','b');

Ts = 1e-5;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','k');

Ts = 0.15
  text(log10(1./Ts),ymax-0.01,'Requirement','rotation',90,...
	'verticalalignment','top','horizontalalignment','right',...
	'fontsize',16,'Parent', h(1));
Ts = 1e-5
  text(log10(1./Ts),ymax-0.01,'E & B','rotation',90,...
	'verticalalignment','top','horizontalalignment','right',...
	'fontsize',16,'fontweight','bold','Parent', h(1));

Ts = 30e-3;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','b','linestyle',':');
text(log10(1./Ts),ymin+0.01,'MMS','rotation',90,...
	'verticalalignment','top','horizontalalignment','left',...
	'fontsize',14,'Parent', h(1));

Ts = 4;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','r','linestyle',':');
text(log10(1./Ts),ymin+0.01,'Cluster','rotation',90,...
	'verticalalignment','top','horizontalalignment','left',...
	'fontsize',14,'Parent', h(1));

Ts = 1/8e3;
plot(h(1),log10(1./Ts).*[1 1],[ymin ymax],'linewidth',linewidth,'color','k','linestyle',':');
text(log10(1./Ts),ymin+0.01,'MMS','rotation',90,...
	'verticalalignment','top','horizontalalignment','left',...
	'fontsize',14,'Parent', h(1));
