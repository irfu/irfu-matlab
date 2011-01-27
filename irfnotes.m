% IRFNOTES File with different common examples how to use irf routines
% use code folding options to fast find your necessary examples
edit irfnotes; return
%% Initializing some figure
% define size to have best agreement with eps file
set(0,'defaultLineLineWidth', 1.5);
fn=figure(61);clf;
set(fn,'color','white'); % white background for figures (default is grey)
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
clear xSize sLeft ySize yTop

% set subplots
% specifying position
h(1)=axes('position',[0.65 0.78 0.2 0.2]); % [x y dx dy]
% having all in standard form
n_subplots=8;i_subplot=1;clear h;
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% 
%% Add information to figures
% text and legends
ht=irf_pl_info([mfilename '  ' datestr(now)]);set(ht,'interpreter','none');

% labels a),b)...
numb={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
for ip=1:2,
  axes(h(ip));
  ht=irf_pl_info(numb{ip},gca,[0.01,1]);
  set(ht,'fontsize',10,'verticalalignment','top');
end
%% Second axis 
hl1 = line(x1,y1,'Color','r');
ax1 = gca;
set(ax1,'XColor','r','YColor','r')

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
%% Reading files
% formatted file reading
%File contents are time intervals in format "T1 T2 Comments":
%2008-03-03T22:50:00 2008-03-03T23:30:00
%2008-03-10T22:10:00 2008-03-10T22:45:00 !
%2008-03-13T07:40:00 2008-03-13T09:40:00 ? shock?

%file reading 
[t1,t2,tint_comments]=textread('Events_reconnection.txt','%s%s%[^\n]');
for j=1:size(t1,1),
  tint(j,1)=iso2epoch(t1{j});tint(j,2)=iso2epoch(t2{j});
end
clear t1 t2 j;
%% Cluster data reading from local disks
% using c_get_batch
    c_get_batch(toepoch([2002 03 04 10 00 00]),30*60,'sp','/home/yuri/caa-data/20020304')
% if time intervals to download are in matrix tint 
for j=1:size(tint,1),
 c_get_batch(tint(j,1),tint(j,2)-tint(j,1),'sp',['./' epoch2iso(tint(j,1),1) '-' epoch2iso(tint(j,2),1)]);
end
clear j;
%% Cluster data reading and plotting from CAA files
for ttt=1, % PANEL: CIS HIA/CODIF spectrogram
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
ic=2;
%dobjname=irf_ssub('C?_CP_CIS_HIA_HS_1D_PEF',ic);eval(['caa_load ' dobjname]);varname=irf_ssub('flux__C?_CP_CIS_HIA_HS_1D_PEF',ic); % HIA
%dobjname=irf_ssub('C?_CP_CIS_CODIF_H1_1D_PEF',ic);eval(['caa_load ' dobjname]);varname=irf_ssub('flux__C?_CP_CIS_CODIF_H1_1D_PEF',ic); % CODIF H+
%dobjname=irf_ssub('C?_CP_CIS_CODIF_O1_1D_PEF',ic);eval(['caa_load ' dobjname]);varname=irf_ssub('flux__C?_CP_CIS_CODIF_O1_1D_PEF',ic); % CODIF O+
%units=eval(['getunits(' dobjname ',''' varname ''')']);
units='log_{10} dEF\newline keV/cm^2 s sr keV';
disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['units: ' units]);
%eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' units ''',''fitcolorbarlabel'');']);
%set(gca,'yscale','log'); set(gca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
%ylabel('E [eV]');
end
for ttt=1, % PANEL: FGM B
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
ic=2;
%dobjname=irf_ssub('C?_CP_FGM_FULL',ic);eval(['caa_load ' dobjname]);varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic);
dobjname=irf_ssub('C?_CP_FGM_5VPS',ic);eval(['caa_load ' dobjname]);varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_5VPS',ic);
c_eval(['B?=getmat(' dobjname ',''' varname ''');'],ic);
c_eval('gsmB?=irf_gse2gsm(B?);',ic);
units=eval(['getunits(' dobjname ',''' varname ''')']);
%units='log_{10} dEF\newline keV/cm^2 s sr keV';
disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['units: ' units]);
c_eval('irf_plot gsmB?',ic);
ylabel('B [nT] GSM');set(gca,'ylim',[-24 19.9]);
irf_legend(gca,{'B_X','B_Y','B_Z'},[0.02 0.35])
%add_timeaxis(gca,'nolabels');
ht=irf_pl_info(['C' num2str(ic)],gca,[0.02,0.7]);set(ht,'fontsize',12,'fontweight','bold')
end

%% CAA Rapid 
dt_rap = 2.0715; % time shift of the pixels

edf = getmat(C1_CP_RAP_ESPCT6,'Electron_Dif_flux__C1_CP_RAP_ESPCT6');
e_e = getmat(C1_CP_RAP_ESPCT6,'Dimension_E__C1_CP_RAP_ESPCT6');

edf(:,2:end) = edf(:,2:end).*e_e'*1.60217646e-12; % convert to erg/cm^2 s sr eV

spec = struct('t',edf(:,1)-dt_rap,'f',e_e(:,1)*1000,'p',[],'f_unit','eV');
spec.p = {edf(:,2:end)};
caa_spectrogram(gca,spec)

set(gca,'YScale','log')
%set(gca,'YTickLabel','$10^5$')
ylabel('')
caxis([-8 -4])
set(h(1),'YLim',[10000.01 245e3],'Color',[.5 .5 .5 ])


% Anisotropy RAPID
e3dd = getmat(C1_CP_RAP_PAD_E3DD,'PAD_Electron_Dif_flux__C1_CP_RAP_PAD_E3DD');

e3dd_data{1} = double(e3dd.data(:,:,1)); % Fix for data gaps in one of the channels
e3dd_data{2} = double(e3dd.data(:,:,5));
e3dd_data{3} = double(e3dd.data(:,:,9));

rap_par = 0.5 * (e3dd_data{1} + e3dd_data{3});
rap_par(e3dd_data{1}==0) = rap_par(e3dd_data{1}==0)*2;
rap_par(e3dd_data{3}==0) = rap_par(e3dd_data{3}==0)*2;

rap_par(rap_par==0) = NaN;

rap_an = e3dd_data{2}./rap_par; rap_an = rap_an - 1;

pcolor(gca,e3dd.t-dt_rap-t_start_epoch,e3dd.dep_x{1}.data(1,:)*1e3,rap_an')
shading flat

caxis(cc);
ylabel('')
set(gca,'TickDir','out','YScale','log','YLim',[10000.01 245e3],'Color',[.5 .5 .5 ])

%% Other
