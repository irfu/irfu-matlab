% IRFNOTES Notes on how to use irf routines
%   Includes different examples that can be directly used
%
%  IRFNOTES   opens file with all the different examples
%
% enable code folding!!! (including cells and if/endif blocks)
% This allows to fast find your necessary examples and execute them.
%

edit irfnotes; return
%% Initialize figure
% fast way
h=irf_plot(5); % h= irf_plot(number_of_subplots);
% example to set subplot position manually
h(1)=axes('position',[0.65 0.78 0.2 0.2]); % [x y dx dy]

% more detailed way
% most lines needed to define the size to have best agreement with eps file
set(0,'defaultLineLineWidth', 1.5);
fn=figure(61);
clf reset;
clear h;
set(fn,'color','white'); % white background for figures (default is grey)
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xSize sLeft ySize yTop
% additional good options
        set(gcf,'defaultAxesFontSize',14);
        set(gcf,'defaultTextFontSize',14);
        set(gcf,'defaultAxesFontUnits','pixels');
        set(gcf,'defaultTextFontUnits','pixels');
%% Print the figure as it looks on screen
set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% to get bitmap file
print -dpng delme.png
% to get pdf file with no white margins
print -depsc2 -painters delme.eps
%print -dpdf -painters delme.pdf
% to convert to pdf on the system command line execute some of
% ps2pdf -dEPSFitPage -dEPSCrop delme.eps
% epstopdf delme.eps
%% Add information to figures
%
% text and legends
help irf_legend
ht=irf_legend(gca,[mfilename '  ' datestr(now)],[0.02 1.01], 'interpreter','none','fontsize',8);
%
% labels a),b)...
help irf_pl_number_subplots
%
% if you want some alternative colorbars, like white in middle, blue negative and red positive
% or if you want to use different colorbars within the same figure
help irf_colormap
%
%% Second axis
%
hl1 = line(x1,y1,'Color','r');
ax1 = gca;
%
ax2 = axes('Position',get(ax1,'Position'));
set(ax2,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
set(ax2,    'YAxisLocation','right');
set(ax2,'Color','none'); % color of axis
set(ax2,'XColor','r','YColor','r'); % color of axis lines and numbers
%


%___________________________________________
% Another example 
% 
t=irf_time([2008 03 01 10 0 0]):.2:irf_time([2008 03 01 11 0 0]);
t=t(:); % make it column vector
% define one time series
y=exp(0.001*(t-t(1))).*sin(2*pi*t/180);
% define data matrix
data1=[t y];
data2=[t sqrt(abs(y))];
% initialize plot
h=irf_plot(1);
% plot data
irf_plot(h(1),data1)
set(h(1),'box','off')
%
h(2) = axes('Position',get(h(1),'Position'));
irf_plot(h(2),data2,'r')
set(h(2),'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
set(h(2),'YAxisLocation','right');
set(h(2),'Color','none','box','off'); % color of axis
set(h(2),'XColor','r','YColor','r'); % color of axis lines and numbers

irf_timeaxis(h(2),'nolabels')

irf_legend(h(1),'data',[0.02 0.98],'color','k')
irf_legend(h(2),'sqrt(data)',[0.8 0.98],'color','r')

ylabel(h(1),'data')
ylabel(h(2),'sqrt(data)')

%% Reading files
% formatted file reading

% txt2mat - excellent routine in matlab exchange central (also put
% into irfu-matlab), see help

%Some example using basic matlab commands
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

%% Example TEMPLATE for a figure
irf_units; % in case you need phyiscal units, see help of irf_units
tint=[irf_time([2006 9 27 17 10 0]) irf_time([2006 9 27 17 40 0])];
ic=1;                % satellite number
CISinstrument='HIA'; % alternative 'CODIF'
if 1, % initialize figure
  fn=figure(61);
  h=irf_plot(8); % number of panels in the plot
  %set(fn,'defaultLineLineWidth',1);  % example of making default changes to all figure
end
if 1, % read the data 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  ADD READ DATA COMMANDS BELOW, see next cells for examples
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
end
if 1, % plot figures panels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  ADD YOUR PANELS BELOW, see next cells for examples
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  
end
if 1, % general commands on all figure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % GENERAL COMMANDS ON ALL FIGURE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  irf_plot_axis_align
  irf_zoom(h,'x',tint);
  irf_legend(h(1),'Figure reference',[0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);
  irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
  irf_timeaxis(h);
end
%% Cluster data reading from CAA, coordinate transformations
%
% To download from CAA
help caa_download
%
% To read downloaded CAA data
tint=[irf_time([2006 9 27 17 14 0]) irf_time([2006 9 27 17 26 0])];
if 0, % read s/c position and velocity
  c_eval('R?=c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');');
  c_eval('V?=c_caa_var_get(''sc_v_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');');
end
if 0, % read FGM data form all sc
 % c_eval('B?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_5VPS'',''mat'');');
  c_eval('B?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');');
  c_eval('B?=irf_abs(B?);');
  c_eval('diB?=c_coord_trans(''GSE'',''ISR2'',B?,''cl_id'',?);');
  c_eval('gsmB?=irf_gse2gsm(B?);');
end
if 0, % read CIS HIA/CODIF velocity moments from available s/c
  c_eval('VCIS?=irf_get_data(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''caa'',''mat'');',[1 3]);
  c_eval('VCISH?=irf_get_data(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''caa'',''mat'');',[1 3 4]);
  c_eval('gsmVCIS?=irf_gse2gsm(VCIS?);',[1 3]);
  c_eval('gsmVCISH?=irf_gse2gsm(VCISH?);',[1 3 4]);
end
if 0, % read RAPID data
  c_eval('[caaRAPID_J?,~,RAPID_J?]=c_caa_var_get(''Electron_Dif_flux__C?_CP_RAP_ESPCT6'');');
end
if 0, % read EFW data
  c_eval('[caaE?,~,diE?]=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'');');
  c_eval('[caaE?,~,diE?]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'');');
  c_eval('[caaVps?,~,Vps?]=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'');');
  c_eval('ExB?=c_caa_var_get(''v_drift_GSE__C?_CP_EFW_L2_V3D_GSE'',''mat'');');
end
if 0, % PEACE calculate density nPEACE1..nPEACE4 [cc] from PITCH_SPIN_DPFlux products above given energy threshold
  for ic=1:4,
    energy_threshold=60; %
    c_eval('caa_load C?_CP_PEA_PITCH_SPIN_DPFlux',ic); % to speed up later
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname);
    phi=c_caa_var_get(var.DEPEND_1);
    x=getfield(c_caa_var_get(phi.DELTA_PLUS),'data');phi_dplus=x(1,:);
    x=getfield(c_caa_var_get(phi.DELTA_MINUS),'data');phi_dminus=x(1,:);
    en=c_caa_var_get(var.DEPEND_2);
    x=getfield(c_caa_var_get(en.DELTA_PLUS),'data');en_dplus=x(1,:);
    x=getfield(c_caa_var_get(en.DELTA_MINUS),'data');en_dminus=x(1,:);
    PEACE_energy_channels=en.data(1,:)+0.5*(en_dplus-en_dminus);
    PEACE_phi_min=phi.data(1,:)-phi_dminus;
    PEACE_phi_max=phi.data(1,:)+phi_dplus;
    ii_energy=find(PEACE_energy_channels>energy_threshold); % use only these energy chanels
    ncoef=ones(length(phi_dplus),length(ii_energy));
    phi_factor=repmat((cos(PEACE_phi_min*pi/180)-cos(PEACE_phi_max*pi/180))',1,length(ii_energy));
    en_factor=repmat(1e-3*(en_dplus(ii_energy)+en_dminus(ii_energy))./sqrt(1e-3*PEACE_energy_channels(ii_energy)),length(phi_dplus),1);
    ncoef=ncoef*0.2284e-7*sqrt(1/1836)*2*pi.*phi_factor.*en_factor;
    
    nPEACE=[varmat.t(:) varmat.t(:)*0];
    varmat.data(isnan(varmat.data))=0;
    for jj=1:size(nPEACE,1),
      nPEACE(jj,2)=sum(sum(shiftdim(varmat.data(jj,:,ii_energy)).*ncoef));
    end
    c_eval('nPEACE?=nPEACE;',ic);
  end
end
if 0, % PEACE calculate density nPEACE1..nPEACE4 [cc] from PITCH_SPIN_DPFlux products using s/c potential correction from EFW
  for ic=1:4,
    c_eval('caa_load C?_CP_PEA_PITCH_SPIN_DPFlux',ic); % to speed up later
    varname=irf_ssub('Spacecraft_potential__C?_CP_EFW_L3_P',ic);
    scpotmat=c_caa_var_get(varname,'mat');
    scpotmat(isnan(scpotmat(:,2)),:)=[]; % remove NaN densities
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    [var,~,varmat]=c_caa_var_get(varname);
    scpot=irf_resamp(scpotmat,varmat); % interpolate sc potential to PEACE data points
    phi=c_caa_var_get(var.DEPEND_1);
    x=getfield(c_caa_var_get(phi.DELTA_PLUS),'data');phi_dplus=x(1,:);
    x=getfield(c_caa_var_get(phi.DELTA_MINUS),'data');phi_dminus=x(1,:);
    en=c_caa_var_get(var.DEPEND_2);
    x=getfield(c_caa_var_get(en.DELTA_PLUS),'data');en_dplus=x(1,:);
    x=getfield(c_caa_var_get(en.DELTA_MINUS),'data');en_dminus=x(1,:);
    PEACE_energy_channels=en.data(1,:)+0.5*(en_dplus-en_dminus);
    PEACE_phi_min=phi.data(1,:)-phi_dminus;
    PEACE_phi_max=phi.data(1,:)+phi_dplus;
    phi_factor=repmat((cos(PEACE_phi_min*pi/180)-cos(PEACE_phi_max*pi/180))',1,length(PEACE_energy_channels));
    
    nPEACE=[varmat.t(:) varmat.t(:)*0];
    varmat.data(isnan(varmat.data))=0;
    for jj=1:size(nPEACE,1),
      satpot=-scpot(jj,2)*1.3; % assumes that probe t spacecraft potential is ~75% of spacecraft potential
      ii_energy=find(PEACE_energy_channels>satpot); % use only these energy chanels
      [en_min,ii_energy_min]=min(PEACE_energy_channels(ii_energy));
      ii_energy(ii_energy==ii_energy_min)=[]; % remove first channel after satellite potential
      [en_min,ii_energy_min]=min(PEACE_energy_channels(ii_energy));
      ii_energy(ii_energy==ii_energy_min)=[]; % remove second channel after satellite potential
      en_factor=repmat(1e-3*(en_dplus(ii_energy)+en_dminus(ii_energy)).*sqrt(1e-3*(PEACE_energy_channels(ii_energy)-satpot))./(1e-3*PEACE_energy_channels(ii_energy)),length(phi_dplus),1);
      ncoef=ones(length(phi_dplus),length(ii_energy));
      ncoef=ncoef*0.2284e-7*sqrt(1/1836)*2*pi.*phi_factor(:,ii_energy).*en_factor;
      nPEACE(jj,2)=sum(sum(shiftdim(varmat.data(jj,:,ii_energy)).*ncoef));
    end
    c_eval('nPEACE?=nPEACE;',ic);
  end
end
if 0, % PEACE calculate pressure P_PEACE1..P_PEACE4 [nPa] from PITCH_SPIN_DPFlux products
  for ic=1:4,
    energy_threshold=60; % only energy channels above this are used (to avoid photoelectrons)
    c_eval('caa_load C?_CP_PEA_PITCH_SPIN_DPFlux',ic); % to speed up later
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname);
    phi=c_caa_var_get(var.DEPEND_1);
    x=getfield(c_caa_var_get(phi.DELTA_PLUS),'data');phi_dplus=x(1,:);
    x=getfield(c_caa_var_get(phi.DELTA_MINUS),'data');phi_dminus=x(1,:);
    en=c_caa_var_get(var.DEPEND_2);
    x=getfield(c_caa_var_get(en.DELTA_PLUS),'data');en_dplus=x(1,:);
    x=getfield(c_caa_var_get(en.DELTA_MINUS),'data');en_dminus=x(1,:);
    PEACE_energy_channels=en.data(1,:)+0.5*(en_dplus-en_dminus);
    PEACE_phi_min=phi.data(1,:)-phi_dminus;
    PEACE_phi_max=phi.data(1,:)+phi_dplus;
    ii_energy=find(PEACE_energy_channels>energy_threshold); % use only these energy chanels
    Pcoef=ones(length(phi_dplus),length(ii_energy));
    phi_factor=repmat((cos(PEACE_phi_min*pi/180)-cos(PEACE_phi_max*pi/180))',1,length(ii_energy));
    en_factor=repmat(1e-3*(en_dplus(ii_energy)+en_dminus(ii_energy)).*sqrt(1e-3*PEACE_energy_channels(ii_energy)),length(phi_dplus),1);
    Pcoef=Pcoef*2/3*0.731026e-8*sqrt(1/1836)*2*pi.*phi_factor.*en_factor;
    
    P_PEACE=[varmat.t(:) varmat.t(:)*0];
    varmat.data(isnan(varmat.data))=0;
    for jj=1:size(P_PEACE,1),
      P_PEACE(jj,2)=sum(sum(shiftdim(varmat.data(jj,:,ii_energy)).*Pcoef));
    end
    c_eval('P_PEACE?=P_PEACE;',ic);
  end
  
end
if 0, % PEACE calculate temperature T_PEACE1..T_PEACE4 [keV] from from nPEACE1..4 and P_PEACE1..4
  for ic=1:4,
    c_eval('T_PEACE?=irf_multiply(4.1609,P_PEACE?,1,nPEACE?,-1);',ic);
  end
end
%
%
%% PANELS that can be used for your figures
% !!! when single sc data are plotted, assumes s/c number is in variable 'ic'

if 1,   % PANEL: C?       FGM Bx,By,Bz,B GSE
  hca=irf_panel('C? FGM B GSE');
  c_eval('irf_plot(hca,B?);',ic);
  ylabel(hca,'B [nT] GSE');
  irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.3])
  irf_legend(hca,{['C' num2str(ic)]},[0.02 0.9],'color','k')
end
if 1,   % PANEL: C?       FGM Bx,By,Bz,B GSM
  hca=irf_panel('C? FGM B GSM');
  c_eval('irf_plot(hca,gsmB?);',ic);
  ylabel(hca,'B [nT] GSM');
  irf_zoom(hca,'y',[-25 25]);
  irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.1])
  irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
if 1,   % PANEL: C?       FGM Bz GSM
  hca=irf_panel('C? FGM Bz');
  c_eval('irf_plot(hca,gsmB?(:,[1 4]));',ic);
  ylabel(hca,'B_Z [nT] GSM');
  irf_legend(hca,{['C' num2str(ic)]},[0.95 0.95],'color','k');
end
if 1,   % PANEL: C1..C4   FGM |B|
  hca=irf_panel(' C1..C4 FGM |B|');
  c_pl_tx(hca,'B?',5)
  ylabel(hca,'|B| [nT]');
  irf_legend(gca,{'C1','C2','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1,   % PANEL: C1..C4   FGM BX GSM
  hca=irf_panel('PANEL: C1..C4, FGM BX');
  c_pl_tx(hca,'gsmB?',2);
  ylabel(hca,'B_X [nT] GSM');
  irf_zoom(hca,'y'); % zoom nicely
  irf_legend(hca,{'C1','C2','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1,   % PANEL: C1..C4   FGM BY GSM
  hca=irf_panel('PANEL: C1..C4, FGM BY');
  c_pl_tx(hca,'gsmB?',3)
  ylabel(hca,'B_Y [nT] GSM');
  irf_legend(hca,{'C1','C2','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1,   % PANEL: C1..C4   FGM BZ GSM
  hca=irf_panel('PANEL: C1..C4, FGM BZ');
  c_pl_tx(hca,'gsmB?',4)
  ylabel(hca,'B_Z [nT] GSM');
  irf_legend(hca,{'C1','C2','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1,   % PANEL: C?       CIS Vx,Vy,Vz,V CODIF(HIA) GSM
  hca=irf_panel('C? CIS V GSM');
  if ic ~=2, % on s/c 2 there is no CIS
    c_eval('irf_plot(hca,gsmVCISH?);',ic);
    % c_eval('irf_plot(hca,gsmVCIS?);',ic); % HIA 
    ylabel(hca,'V [km/s] GSM');
    irf_zoom(hca,'y',[-200 1200]);
    irf_legend(hca,{'V_X','V_Y','V_Z'},[0.02 0.49])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
  end
end
if 1,   % PANEL: C1,C3,C4 CIS Vx GSM velocities
  hca=irf_panel(h,'C1,C3,C4 CIS Vx velocities');
  hold(hca,'off');
  irf_plot(hca,gsmVCIS1(:,1:2),'color','k'); % HIA
  hold(hca,'on');
  irf_plot(hca,gsmVCIS3(:,1:2),'color','g'); % HIA
  irf_plot(hca,gsmVCISH4(:,1:2),'color','b'); % CODIF
  ylabel(hca,'V_X [km/s] GSM');
  irf_zoom(hca,'y',[-1000 1000]);
  irf_legend(hca,{'C1','','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1,   % PANEL: C1,C3,C4 CIS Vy GSM velocities
  hca=irf_panel(h,'C1,C3,C4 CIS Vy velocities');
  hold(hca,'off');
  irf_plot(hca,gsmVCIS1(:,[1 3]),'color','k'); % HIA
  hold(hca,'on');
  irf_plot(hca,gsmVCIS3(:,[1 3]),'color','g'); % HIA
  irf_plot(hca,gsmVCISH4(:,[1 3]),'color','b'); % CODIF
  ylabel(hca,'V_X [km/s] GSM');
  irf_zoom(hca,'y',[-1000 1000]);
  irf_legend(hca,{'C1','','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1,   % PANEL: C1,C3,C4 CIS Vz GSM velocities
  hca=irf_panel(h,'C1,C3,C4 CIS Vz velocities');
  hold(hca,'off');
  irf_plot(hca,gsmVCIS1(:,[1 4]),'color','k'); % HIA
  hold(hca,'on');
  irf_plot(hca,gsmVCIS3(:,[1 4]),'color','g'); % HIA
  irf_plot(hca,gsmVCISH4(:,[1 4]),'color','b'); % CODIF
  ylabel(hca,'V_Z [km/s] GSM');
  irf_zoom(hca,'y',[-150 150]);
  irf_legend(hca,{'C1','','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1,   % PANEL: Pressures, B and CIS HIA/CODIF single s/c
  hca=irf_subplot('C? CIS/FGM Pressures')
  if ic~=2,
    if strcmp(CISinstrument,'HIA')
      dobjname=irf_ssub('C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic);
      caa_load(dobjname);
      varname=irf_ssub('pressure__C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic);
      c_eval(['PressureCIS?=getmat(' dobjname ',''' varname ''');'],ic); % Pressure in nPa
      varname=irf_ssub('density__C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic);
      c_eval(['DensityCIS?=getmat(' dobjname ',''' varname ''');'],ic); % Density in cc
      varname=irf_ssub('velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic);
      c_eval(['VelocityCIS?=getmat(' dobjname ',''' varname ''');'],ic); % Velocity in km/s
    elseif strcmp(CISinstrument,'CODIF')
      dobjname=irf_ssub('C?_CP_CIS_CODIF_HS_H1_MOMENTS',ic);
      caa_load(dobjname);
      varname=irf_ssub('velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS',ic);
      c_eval(['VelocityCIS?=getmat(' dobjname ',''' varname ''');'],ic); % Velocity in km/s
      varname=irf_ssub('density__C?_CP_CIS_CODIF_HS_H1_MOMENTS',ic);
      c_eval(['DensityCIS?=getmat(' dobjname ',''' varname ''');'],ic); %
      varname=irf_ssub('T__C?_CP_CIS_CODIF_HS_H1_MOMENTS',ic);
      c_eval(['TemperatureCIS?=getmat(' dobjname ',''' varname ''');'],ic); %
      c_eval('PressureCIS?=irf_multiply(Units.kB*1e6*1e6*1e9,DensityCIS?,1,TemperatureCIS?,1);',ic); %
    end
    c_eval('VelocityCIS?=irf_abs(VelocityCIS?);',ic); %
    c_eval('DynamicPressureCIS?=irf_multiply(0.5*1e6*Units.mp*1e6*1e9,DensityCIS?,1,VelocityCIS?(:,[1 5]),2);',ic); % in nPa, assumes H+
    c_eval('PressureB?=irf_multiply(1e-18*0.5/Units.mu0*1e9,B?(:,[1 5]),2,B?(:,[1 5]),0);',ic); % Pressure in nPa
    c_eval('PressureTotal?=irf_add(1,PressureCIS?,1,PressureB?);',ic); % Pressure in nPa
    varunits=eval(['getunits(' dobjname ',''' varname ''')']);
    %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
    disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
    c_eval('irf_plot(hca,{PressureB?,PressureCIS?,DynamicPressureCIS?,PressureTotal?},''comp'');',ic);
    ylabel(hca,'P [nPa]');
    set(hca,'ylim',[0 0.29]);
    irf_legend(hca,{'P_B','P_i','P_{kin}'},[0.02 0.3])
    irf_legend(hca,{['C' num2str(ic)]},[0.02 0.9],'color','k')
  end
end
irf_colormap(hca,'default'); % execute once to change the colormap or for every axis you want to fix the different colormap
if 1,   % PANEL: C?       CIS HIA/CODIF spectrogram
  hca=irf_panel('C? CIS HIA/CODIF spectrogram');
  if ic~=2,
    if ic==4 || strcmp(CISinstrument,'CODIF')
      varname=irf_ssub('flux__C?_CP_CIS_CODIF_H1_1D_PEF',ic); % CODIF H+
    elseif strcmp(CISinstrument,'HIA')
      varname=irf_ssub('flux__C?_CP_CIS_HIA_HS_1D_PEF',ic); % HIA
    end
    %varunits=eval(['getunits(' dobjname ',''' varname ''')']);
    varunits={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
    caxis(hca,[3.9 6.1]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E [eV]');
  end
  irf_legend(hca,{['C' num2str(ic)]},[0.98 0.1],'color','k')
end
if 1,   % PANEL: C?       EFW Ex,Ey ISR2 
  hca=irf_panel('EFW E ISR2 c?');
  c_eval('irf_plot(hca,diE?(:,1:3))',ic);
  ylabel(hca,'E [mV/m] ISR2');
  irf_zoom(hca,'ylim','smart');
  irf_legend(hca,{'E_X','E_Y'},[0.02 0.98])
  irf_legend(hca,{['C' num2str(ic)]},[0.02 0.95],'color','k')
end
if 1,   % PANEL: C?       EFW satellite potential
  % change L3 to L2 to get full resolution instead of spin
  hca=irf_panel('EFW satellite potential spin');
  Vps=irf_ssub('Spacecraft_potential__C?_CP_EFW_L3_P',ic);
  irf_plot(hca,Vps,'nolabels');
  irf_zoom(hca,'y',[-50 -10]);
  ylabel(hca,'Sat pot [V]');
end
if 1,   % PANEL: C?       EFW ExB - velocity
  % change L3 to L2 to get full resolution instead of spin
  hca=irf_panel('EFW satellite potential spin');
  ExB=irf_ssub('v_drift_GSE__C?_CP_EFW_L3_V3D_GSE',ic);
  irf_plot(hca,ExB,'nolabels');
  irf_zoom(hca,'y',[-500 500]);
  ylabel(hca,'ExB [km/s]');
end
if 0,   % PANEL: C1..C4   EFW satellite potential
  % change L3 to L2 to get full resolution instead of spin
  hca=irf_panel('EFW satellite potential spin');
  hold(hca,'off');
  irf_plot(hca,'Spacecraft_potential__C1_CP_EFW_L3_P');
  hold(hca,'on');
  irf_plot(hca,'Spacecraft_potential__C2_CP_EFW_L3_P','color','r');
  irf_plot(hca,'Spacecraft_potential__C3_CP_EFW_L3_P','color',[0 0.5 0]);
  irf_plot(hca,'Spacecraft_potential__C4_CP_EFW_L3_P','color','b');
  
  irf_zoom(hca,'y',[-50 -10]);
  irf_legend(hca,{'C1','C2','C3','C4'},[0.98, 0.1],'color','cluster');
  ylabel(hca,'Sat pot [V]');
end
if 1,   % PANEL: C?       STAFF spectrogram Bx
  hca=irf_panel('STAFF spectrogram B and fce/flh lines');
  varname=irf_ssub('BB_xxyyzz_isr2__C?_CP_STA_PSD',ic);
  [~,dobj,~,varunits]=c_caa_var_get(varname);
  %varunits='nT^2/Hz';
  hold(hca,'off');
  plot(hca,dobj,varname,'colorbarlabel',varunits,'fitcolorbarlabel','comp',1);
  hold(hca,'on');
  c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic); % calculate electron gyrofrequency
  irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
  irf_plot(hca,[fce(:,1) fce(:,2)*.5],'-','linewidth',0.2,'color','w'); % fce/2 white line
  irf_plot(hca,[fce(:,1) fce(:,2)*.25],'-','linewidth',0.2,'color','w'); % fce/4 white line
  c_eval('flh=irf_plasma_calc(B?,1,0,0,0,''Flh'');',ic); % calculate lower hybrid frequency (underdense case in puter magnetosphere)
  irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
  caxis(hca,[-10 -7]);
  set(hca,'yscale','log','ytick',[1e1 1e2 1e3]);
  irf_zoom(hca,'y',[10 4000]);
end
if 1,   % PANEL: C?       STAFF spectrogram Ex and fce/flh lines
  hca=irf_panel('STAFF spectrogram Ex and fce/flh lines');
  varname=irf_ssub('EE_xxyy_isr2__C?_CP_STA_PSD',ic);
  %[~,~,~,varunits]=c_caa_var_get(varname);
  varunits='(mV/m)^2/Hz';
  irf_plot(hca,varname,'tint',tint,'colorbarlabel',varunits,'fitcolorbarlabel','comp',1,'nolabels');
  % next lines are examples how to add gyro/lower hybrid frequency lines
  % B & n should be calculated before
  hold(hca,'on');
  c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic); % calculate electron gyrofrequency
  irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
  c_eval('fpe=irf_plasma_calc(irf_resamp(B?,ncal_PEACE?),ncal_PEACE?,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
  irf_plot(hca,fpe,'-','linewidth',0.2,'color','k');
  c_eval('fpemom=irf_plasma_calc(irf_resamp(B?,nPEACE),nPEACE,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
  irf_plot(hca,fpemom,'.','linewidth',0.2,'color','r');
  c_eval('flh=irf_plasma_calc(B?,1,0,0,0,''Flh'');',ic); % calculate lower hybrid frequency (underdense case in puter magnetosphere)
  irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
  hold(hca,'off');
  % polish the panel
  caxis(hca,[-9 -1]);
  set(hca,'yscale','log','ytick',[1e1 1e2 1e3]);
  irf_zoom(hca,'y',[10 4000]);
end
if 1,   % PANEL: C?       WHISPER spectrogram
  hca=irf_panel('WHISPER spectrogram natural');
  varname=irf_ssub('Electric_Spectral_Power_Density__C?_CP_WHI_NATURAL',ic);
  varunits='(V/m)^2/Hz';
  % REMOVE 'fillspectgrogramgaps' flag in the next line if precise intervals of
  % WHISPER measurements are needed !!!!
  % If working with shorter intervals can also remove 'tint' option
  irf_plot(hca,varname,'tint',tint,'colorbarlabel',varunits,'fitcolorbarlabel','fillspectrogramgaps','nolabels');
  if 1, % add plasma frequency lines
      % density n should be calculated before
      hold(hca,'on');
      c_eval('fpe=irf_plasma_calc(irf_resamp(1,n),n,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
      irf_plot(hca,irf_tappl(fpe,'/1e3'),'-','linewidth',0.2,'color','k');
      hold(hca,'off');
  end
  % polish the panel
  caxis(hca,[-16 -11]);
  set(hca,'yscale','log','ytick',[3 4 5 1e1 20 30 50 ]);
  irf_zoom(hca,'y',[2 12]);
end
if 1,   % PANEL: C?       RAPID electron spectrogram
  hca=irf_panel('C?_CP_RAP_ESPCT6');
  varname=irf_ssub('Electron_Dif_flux__C?_CP_RAP_ESPCT6',ic);
  %varunits=irf_get_data(varname,'caa','unit');
  varunits={'log_{10} dF','1/cm^2 s sr keV'};
  irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
  caxis(hca,[0.51 4.49]);
  ylabel(hca,'E [keV]');
  set(hca,'yscale','log');
  set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
  irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
if 1,   % PANEL: C?       RAPID ion spectrogram
    hca=irf_panel('C?_CP_RAP_HSPCT');
    varname=irf_ssub('Proton_Dif_flux__C?_CP_RAP_HSPCT',ic);
    %varunits=c_caa_var_get(varname,'units');
    varunits={'log_{10} dF','1/cm^2 s sr keV'};
    irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
    caxis(hca,[0.51 4.49]);
    set(hca,'yscale','log');
    ylabel(hca,'E [keV]');
    irf_zoom(hca,'y',[80 2000]);
    set(hca,'ytick',[1e2 2e2 5e2 1e3 1e4 1e5])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
if 0,   % PANEL: RAPID spectrogram parallel
  hca=h(i_subplot);i_subplot=i_subplot+1;
  dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
  caa_load(dobjname);
  varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
  varunits=eval(['getunits(' dobjname ',''' varname ''')']);
  l3dd_caa=eval(['getv(' dobjname ',''' varname ''')']);
  l3dd=eval(['getmat(' dobjname ',''' varname ''')']);
  l3dd_energies=eval(['getv(' dobjname ',''' l3dd_caa.DEPEND_1 ''')']);
  disp(['PANEL: C' num2str(ic) ' RAPID anisotropy']);
  disp(['dobj:' dobjname ]);disp([' var:' varname]);
  eval(['plot(hca,' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',1);']);
  caxis(hca,[1.1 4.3]);
  irf_colormap;
  set(hca,'yscale','log');
  set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5]);
end
if 0,   % PANEL: RAPID PAD_L3DD spectrogram perpendicular
  hca=h(i_subplot);i_subplot=i_subplot+1;
  dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
  caa_load(dobjname);
  varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
  l3dd_caa=eval(['getv(' dobjname ',''' varname ''')']);
  l3dd=eval(['getmat(' dobjname ',''' varname ''')']);
  l3dd_energies=eval(['getv(' dobjname ',''' l3dd_caa.DEPEND_1 ''')']);
  disp(['PANEL: C' num2str(ic) ' RAPID anisotropy']);
  disp(['dobj:' dobjname ]);disp([' var:' varname]);
  eval(['plot(hca,' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',5);']);
  caxis(hca,[1.1 4.3]);
  irf_colormap;
  set(hca,'yscale','log');
  set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5]);
end
if 0,   % PANEL: RAPID spectrogram pitch angle
  hca=h(i_subplot);i_subplot=i_subplot+1;
  dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
  caa_load(dobjname);
  varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
  disp(['PANEL: C' num2str(ic) ' RAPID anisotropy']);
  disp(['dobj:' dobjname ]);disp([' var:' varname]);
  eval(['plot(hca,' dobjname ',''' varname ''',''ax'',hca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp_dim1'',''comp'',1);']);
  caxis(hca,[1.1 4.3]);
  irf_colormap;
  set(hca,'yscale','lin');
  set(hca,'ytick',[0 45 90 135 180]);
  ylabel(hca,'Pitch ang. [deg]');
end
if 0,   % PANEL: RAPID spectrogram anisotropy
  hca=h(i_subplot);i_subplot=i_subplot+1;
  dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
  caa_load(dobjname);
  varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
  l3dd_caa=eval(['getv(' dobjname ',''' varname ''')']);
  l3dd=eval(['getmat(' dobjname ',''' varname ''')']);
  l3dd_energies=eval(['getv(' dobjname ',''' l3dd_caa.DEPEND_1 ''')']);
  disp(['PANEL: C' num2str(ic) ' RAPID anisotropy']);
  disp(['dobj:' dobjname ]);disp([' var:' varname]);
  
  % Fix for data gaps in one of the channels
  rap_0deg = double(l3dd.data(:,:,1));
  rap_perp = double(l3dd.data(:,:,5));
  rap_180deg = double(l3dd.data(:,:,9));
  
  rap_par = 0.5 * (rap_0deg + rap_180deg);
  rap_par(rap_0deg==0) = rap_par(rap_0deg==0)*2;
  rap_par(rap_180deg==0) = rap_par(rap_180deg==0)*2;
  
  rap_par(rap_par==0) = NaN;
  
  rap_an = rap_perp./rap_par; rap_an = rap_an - 1;
  % check if DELTA_PLUS and  DELTA_MINUS are given
  if isfield(l3dd_energies,'DELTA_PLUS') && isfield(l3dd_energies,'DELTA_MINUS')
    specrec.df=struct('plus',l3dd_energies.DELTA_PLUS,'minus',l3dd_energies.DELTA_MINUS);
    if ischar(l3dd_energies.DELTA_PLUS)
      deltaplus= getv(dobj,l3dd_energies.DELTA_PLUS);
      deltaminus= getv(dobj,l3dd_energies.DELTA_MINUS);
      dep_x{d}.df.plus=deltaplus.data(1,:);
      dep_x{d}.df.minus=deltaminus.data(1,:);
    end
  else specrec.df=[];
  end
  
  specrec.t=l3dd.t;
  specrec.f=l3dd_energies.data(1,:);
  specrec.p={rap_an};
  irf_spectrogram(hca,specrec);
  colorbar('peer',hca);
  caxis(hca,[-2 2]);
  irf_colormap('poynting');
  set(hca,'yscale','log');
  set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
end
if 0,   % PANEL: C1..C4   PEACE density MOMENTS
  hca=irf_panel('PEACE N MOMENTS');
  c_eval('nPEACE?=irf_get_data(''Data_Density__C?_CP_PEA_MOMENTS'',''caa'',''mat'');');
  c_pl_tx(hca,'nPEACE?','-');
  set(hca,'yscale','log','ytick',[1e-2 1e-1 1e0 1e1]);
  irf_zoom(hca,'y',[0.005 0.3]);
  irf_legend(hca,{'C1','C2','C3','C4'},[0.98, 0.1],'color','cluster');
  ylabel(hca,'N_{e} [cc]');
end
if 0,   % PANEL: C1..C4   PEACE density MOMENTS
  hca=irf_panel('C1..C4 PEACE T MOMENTS');
  c_eval('T_PEACE?=irf_get_data(''Data_Temperature_ComponentPerpendicularToMagField__C2_CP_PEA_'',''caa'',''mat'');');
  c_eval('T_PEACE?=irf_tappl(T_PEACE?,''*86.132'');'); % conversion from MK to eV 
  c_pl_tx(hca,'T_PEACE?','.-');
  set(hca,'yscale','log','ytick',[1e-2 1e-1 1e0 1e1]);
%  irf_zoom(hca,'y',[0.005 0.3]);
  irf_legend(hca,{'C1','C2','C3','C4'},[0.98, 0.1],'color','cluster');
  ylabel(hca,'T_{e} [eV]');
end
if 0,   % PANEL: C?       PEACE PITCH_SPIN_DEFlux spectrogram omni
    hca=irf_panel('C? PEACE energy spectra');
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
    %varunits=irf_get_data(varname,'caa','units');
    varunits={'log_{10} dEF','keV/cm^2 s sr keV'};
    irf_plot(hca,varname,'sum_dim1','colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
    caxis(hca,[5.6 7.9]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E [eV]');
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.2],'color','k')
end
if 0,   % PANEL: C?       PEACE PITCH_SPIN_DEFlux spectrogram angles
    hca=irf_panel('C? PEACE DEFlux pitch spectra');
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
    varunits=irf_get_data(varname,'caa','units');
    %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
    irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','comp',22);
    caxis(hca,[5.9 8.9]);
    set(hca,'yscale','lin');
    set(hca,'ytick',[ 45 90 135 ]);
    ylabel(hca,'\Theta [deg]');
end
if 0,   % PANEL: C?       PEACE PITCH_SPIN_PSD spectrogram angles
    hca=irf_panel('C? PEACE PSD pitch spectra');
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_PSD',ic);
    varunits=irf_get_data(varname,'caa','units');
    %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
    irf_plot(hca,varname,'colorbarlabel',varunits,'fitcolorbarlabel','comp',22);
    caxis(hca,[1.4 3.1]);
    set(hca,'yscale','lin');
    set(hca,'ytick',[ 45 90 135 ]);
    ylabel(hca,'\Theta [deg]');
end
if 0,   % PANEL: C?       PEACE PITCH_SPIN_DEFlux spectrogram angles in specified energy range
    hca=irf_panel('C? PEACE DEFlux pitch spectra 2');
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
    pea=c_caa_var_get(varname,'mat');
    Enmin=1e3;Enmax=5e3; % energy interval
    energy_channels=[]; % if do not want energy interval, specify energy channel
    energy_channel_values=pea.dep_x{2}.data(1,:);
    if numel(energy_channels)==0, % specified energy interval
        energy_channels=find(energy_channel_values>Enmin & energy_channel_values < Enmax); % use only these energy chanels
        en_values=energy_channel_values(energy_channels);
        en_label=[num2str(min(en_values),'%5.0f') ' - ' num2str(max(en_values),'%5.0f') ' ' pea.dep_x{2}.units];
    elseif numel(energy_channels)>1,
        en_values=energy_channel_values(energy_channels);
        en_label=[num2str(min(en_values),'%5.0f') ' - ' num2str(max(en_values),'%5.0f') ' ' pea.dep_x{2}.units];
    elseif numel(energy_channels)==1, % only one channel specified
        en_label=[num2str(energy_channel_values(energy_channels)) ' ' pea.dep_x{2}.units];
    end
    specrec={};
    specrec.p=squeeze(sum(pea.data(:,:,energy_channels),3));
    specrec.t=pea.t;
    specrec.f=pea.dep_x{1}.data(1,:);
    specrec.dt=pea.dt;
    specrec.df=pea.dep_x{1}.df;
    varunits=irf_get_data(varname,'caa','units');
    specrec.p_label=varunits;
    %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
    irf_spectrogram(hca,specrec);
    irf_legend(hca,en_label,[0.98 0.05]);
    irf_legend(hca,['C' num2str(ic)],[0.98 0.98]);
    caxis(hca,[5.9 8.9]);
    set(hca,'yscale','lin');
    set(hca,'ytick',[ 45 90 135 ]);
    ylabel(hca,'\Theta [deg]');
end
if 0,   % PANEL: RAPID/PEACE electron densities
  % requires that both RAPID and PEACE densities are calculated before
  hca=irf_panel('nonthermal_vs_thermal_e');
  %
  % RAPID density
  %
  varname=irf_ssub('Electron_Dif_flux__C?_CP_RAP_ESPCT6',ic);
  [var,~,varmat,varunits]=c_caa_var_get(varname);
  en=c_caa_var_get(var.DEPEND_1);
  en_dplus=en.DELTA_PLUS;
  RAPID_energy_channels=en.data(1,:)+0.5*en_dplus;
  ncoef=0.2284e-7*sqrt(1/1836)*4*pi*en_dplus./sqrt(RAPID_energy_channels);
  nmat=varmat; nmat(:,2:end)=nmat(:,2:end).*repmat(ncoef,size(varmat,1),1);
  nRAPID=nmat(:,1:2);nRAPID(:,2)=sum(nmat(:,2:end),2);
  %
  % PEACE density
  %
  varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
  caa_load('C?_CP_PEA_PITCH_SPIN_DPFlux'); % speeds up later loading
  [var,dobj,varmat,varunits]=c_caa_var_get(varname);
  phi=c_caa_var_get(var.DEPEND_1);
  x=getfield(c_caa_var_get(phi.DELTA_PLUS),'data');phi_dplus=x(1,:);
  x=getfield(c_caa_var_get(phi.DELTA_MINUS),'data');phi_dminus=x(1,:);
  en=c_caa_var_get(var.DEPEND_2);
  x=getfield(c_caa_var_get(en.DELTA_PLUS),'data');en_dplus=x(1,:);
  x=getfield(c_caa_var_get(en.DELTA_MINUS),'data');en_dminus=x(1,:);
  PEACE_energy_channels=en.data(1,:)+0.5*(en_dplus-en_dminus);
  PEACE_phi_min=phi.data(1,:)-phi_dminus;
  PEACE_phi_max=phi.data(1,:)+phi_dplus;
  ii_energy=find(PEACE_energy_channels>100); % use only these energy chanels
  ncoef=ones(length(phi_dplus),length(ii_energy));
  phi_factor=repmat((cos(PEACE_phi_min*pi/180)-cos(PEACE_phi_max*pi/180))',1,length(ii_energy));
  en_factor=repmat(1e-3*(en_dplus(ii_energy)+en_dminus(ii_energy))./sqrt(1e-3*PEACE_energy_channels(ii_energy)),length(phi_dplus),1);
  ncoef=ncoef*0.2284e-7*sqrt(1/1836)*2*pi.*phi_factor.*en_factor;
  nPEACE=[varmat.t(:) varmat.t(:)*0];
  varmat.data(isnan(varmat.data))=0;
  for jj=1:size(nPEACE,1),
    nPEACE(jj,2)=sum(sum(shiftdim(varmat.data(jj,:,ii_energy)).*ncoef));
  end
  %
  % Ratio of densities
  %
  ratio_RAPIDvsPEACE=irf_multiply(1e3,nRAPID,1,nPEACE,-1);
  %
  % Plot
  %
  irf_plot(hca,ratio_RAPIDvsPEACE,'.-');
  ylabel(hca,'(n_{RAPID}/n_{PEACE})10^3')
  set(hca,'yscale','lin');
end
if 0,   % PANEL: RAPID PSD power law fit
  hca=h(i_subplot); i_subplot=i_subplot+1;
  [caaFlux_RAP,~,Flux_RAP]=c_caa_var_get('PAD_Electron_Dif_flux__C2_CP_RAP_PAD_E3DD');
  DPF2PSD=[5.369 4.487 3.872 3.517 3.402 3.556 4.080 4.883].*1e-9;
  PSD_RAP=Flux_RAP;
  for ii=1:8,
    PSD_RAP.data(:,ii,:)=Flux_RAP.data(:,ii,:) * DPF2PSD(ii);
  end
  nspins_to_average=2; % how many spins to average for fit
  n_pitchangles_to_average_par=[1 2 8 9]; % pitch angles to average , parallel psd
  n_pitchangles_to_average_perp=[4 5 6];  % pitch angles to average for perp psd
  for jj=nspins_to_average:nspins_to_average:size(PSD_RAP.data,1),
    ind=jj-nspins_to_average+1:jj;
    xx=shiftdim(sum(PSD_RAP.data(ind,:,:),1)/nspins_to_average); % time average
    xx_par=sum(xx(:,n_pitchangles_to_average_par),2)/numel(n_pitchangles_to_average_par);
    xx_perp=sum(xx(:,n_pitchangles_to_average_perp),2)/numel(n_pitchangles_to_average_perp);
    PSD_RAP.psdpar(jj/nspins_to_average,:)=xx_par;
    PSD_RAP.psdperp(jj/nspins_to_average,:)=xx_perp;
    PSD_RAP.tav(jj/nspins_to_average)=sum(PSD_RAP.t(ind))/nspins_to_average;
  end
  PSD_RAP.tav=PSD_RAP.tav(:); % to get column vector
  PSD_RAP.psdpar=log10(PSD_RAP.psdpar);
  PSD_RAP.psdperp=log10(PSD_RAP.psdperp);
  PSD_RAP.kpar=zeros(size(PSD_RAP.psdpar,1),1);   % allocate matrix
  PSD_RAP.kperp=zeros(size(PSD_RAP.psdperp,1),1); % allocate matrix
  en=c_caa_var_get('Dimension_E__C2_CP_RAP_PAD_E3DD');
  energy_levels=en.data(1,:)+0.5*en.DELTA_PLUS;
  log10_energy_levels=log10(energy_levels);
  for jj=1:size(PSD_RAP.psdpar,1),
    ind_noninf=~isinf(PSD_RAP.psdpar(jj,:)); % use only point with counts (log(zero counts)=-Inf)
    [p,s]=polyfit(log10_energy_levels(ind_noninf),PSD_RAP.psdpar(jj,ind_noninf),1);
    PSD_RAP.kpar(jj)=p(1);
    ind_noninf=~isinf(PSD_RAP.psdperp(jj,:)); % use only point with counts (log(zero counts)=-Inf)
    [p,s]=polyfit(log10_energy_levels(ind_noninf),PSD_RAP.psdperp(jj,ind_noninf),1);
    PSD_RAP.kperp(jj)=p(1);
  end
  irf_plot(hca,[PSD_RAP.tav PSD_RAP.kpar PSD_RAP.kperp ]);
  irf_zoom([-6 -2],'y',hca);
  irf_legend(hca,{'k_{par}','k_{perp}'},[0.95 0.85]);
  ylabel(hca,'Power law slope');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subspin resolution panels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0,   % PANEL: PEACE PEA_PITCH_3DRH_PSD high res
  hca=irf_panel('C? PEACE 3DRH');
  ic=3;
  res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DRH_PSD',ic));
  [delmett,ind]=irf_tlim(res.tt,tint);
  specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
  if 0, % energy spectorgram (integrated over pitch angles)
    specrec.f=log10(res.en);
    specrec.p=res.omni(ind,:);
    specrec.f_label=['Log10 ' res.enlabel];
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    enindex=15;
    res.en(enindex)
    specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
    specrec.p=log10(res.data(ind,:,enindex));
  end
  specrec_C33DRH=specrec;
  irf_spectrogram(hca,specrec);
  caxis(hca,[-2.99 0.99])
  irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
  set(hca,'ytick',[30 60 90 120 150]);
end
if 1,   % PANEL: PEACE 3DXPH_DEFlux high res energy spectrogram
  hca=irf_panel('PEACE 3DXPH_DEFlux energy');
  res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_3DXPH_PSD',ic));
  [~,ind]=irf_tlim(res.tt,tint);
  specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
  if 1, % energy spectorgram (integrated over pitch angles)
    specrec.f=log10(res.en);
    specrec.p=res.omni(ind,:);
    specrec.f_label=['Log10 ' res.enlabel];
    irf_spectrogram(hca,specrec);
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    enindex=13;
    specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
    specrec.p=log10(res.data(ind,:,enindex));
    irf_spectrogram(hca,specrec);
    set(hca,'ytick',[30 60 90 120 150]);
  end
  caxis(hca,[-1.99 0.49]);
  irf_legend(hca,['C' num2str(ic)],[0.98,0.98]);
end
if 1,   % PANEL: PEACE 3DXPH_DEFlux high res angular spectrogra,
  hca=irf_panel('PEACE 3DXPH_DEFlux angular');
  res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_3DXPH_PSD',ic));
  [delmett,ind]=irf_tlim(res.tt,tint);
  specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
  if 0, % energy spectorgram (integrated over pitch angles)
    specrec.f=log10(res.en);
    specrec.p=res.omni(ind,:);
    specrec.f_label=['Log10 ' res.enlabel];
    irf_spectrogram(hca,specrec);
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    enindex=13; % specify which energy chanel
    specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
    specrec.p=log10(res.data(ind,:,enindex));
    irf_spectrogram(hca,specrec);
    set(hca,'ytick',[30 60 90 120 150]);
  end
  caxis(hca,[-1.99 0.49]);
  irf_legend(hca,['C' num2str(ic)],[0.98,0.98]);
end
if 0,   % PANEL: RAPID L3DD high res pitch C4
  h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);hca=h(i_subplot);i_subplot=i_subplot+1;
  ic=4;
  res=c_caa_construct_subspin_res_data(irf_ssub('Electron_L_Dif_flux_3D__C?_CP_RAP_L3DD',ic));
  [delmett,ind]=irf_tlim(res.tt,tint);
  specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label','Log PSD [s^3/km^6]');
  if 0, % energy spectrogram (integrated over pitch angles)
    specrec.f=res.en;
    specrec.p=res.omni(ind,:);
    specrec.f_label=[res.enlabel];
    yticks=[1 2 3 4 5];
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    enindex=1;specrec.f_label=[specrec.f_label '\newline  [E=' num2str(res.en(enindex),4) 'keV]'];
    specrec.p=log10(res.data(ind,:,enindex)*5.369e-9);
    yticks=[30 60 90 120 150];
  end
  irf_spectrogram(hca,specrec);
  caxis(hca,[-4.99 -3.01])
  irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
  set(hca,'ytick',yticks);
end
if 0,   % PANEL: CIS CODIF high res energy C4
  h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);hca=h(i_subplot);i_subplot=i_subplot+1;
  res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_CODIF_HS_H1_PSD',ic));
  %res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_HIA_HS_MAG_IONS_PSD',ic));
  
  [delmett,ind]=irf_tlim(res.tt,tint);
  specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
  if 0, % energy spectorgram (integrated over pitch angles)
    specrec.f=log10(res.en);
    specrec.p=res.omni(ind,:);
    specrec.f_label='Log_{10} E [eV]';
    yticks=[1 2 3 4 5];
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    enindex=(26:30);
    if numel(enindex)==1,
      specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex),'%6.f') 'eV]'];
    else
      specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex(1)),'%6.f') ' - ' num2str(res.en(enindex(end)),'%6.f') 'eV]'];
    end
    specrec.p=sum(res.data(ind,:,enindex),3);
    yticks=[45 90 135 ];
  end
  irf_spectrogram(hca,specrec);colormap(jet);
  caxis([1 5]);
  set(hca,'ytick',yticks);
  irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
end
if 1,   % PANEL: CIS HIA/CODIF high res energy C3
  hca=irf_panel('CIS CODIF high res energy');
  ic=3;
  %res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_CODIF_HS_H1_PSD',ic));
  res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_HIA_HS_MAG_IONS_PEF',ic));
  
  [~,ind]=irf_tlim(res.tt,tint);
  specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PEF \newline [' res.dataunits ']']);
  if 1, % energy spectorgram (integrated over pitch angles)
    specrec.f=log10(res.en);
    specrec.p=res.omni(ind,:);
    specrec.f_label='Log_{10} E [eV]';
    yticks=[1 2 3 4 5];
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    enindex=(26:30);
    if numel(enindex)==1,
      specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex),'%6.f') 'eV]'];
    else
      specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex(1)),'%6.f') ' - ' num2str(res.en(enindex(end)),'%6.f') 'eV]'];
    end
    specrec.p=sum(res.data(ind,:,enindex),3);
    yticks=[45 90 135 ];
  end
  irf_spectrogram(hca,specrec);
  caxis([1 5]);
  set(hca,'ytick',yticks);
  irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
end
if 1,   % PANEL: CIS HIA/CODIF high res energy C4
  hca=irf_panel('CIS CODIF high res pitch angle');
  %res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_CODIF_HS_H1_PSD',ic));
  res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_HIA_HS_MAG_IONS_PEF',ic));
  
  [~,ind]=irf_tlim(res.tt,tint);
  specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PEF  \newline [' res.dataunits ']']);
  if 0, % energy spectorgram (integrated over pitch angles)
    specrec.f=log10(res.en);
    specrec.p=res.omni(ind,:);
    specrec.f_label='Log_{10} E [eV]';
    yticks=[1 2 3 4 5];
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=res.theta;specrec.f_label='Pitch angle';
    specrec.p=res.pitch_angle(ind,:);
    enindex=(26:30);
    if numel(enindex)==1,
      specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex),'%6.f') 'eV]'];
    else
      specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex(1)),'%6.f') ' - ' num2str(res.en(enindex(end)),'%6.f') 'eV]'];
    end
    specrec.p=sum(res.data(ind,:,enindex),3);
    yticks=[45 90 135 ];
  end
  irf_spectrogram(hca,specrec);
  caxis([4 7]);
  set(hca,'ytick',yticks);
  irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
end

%% Another example plot, simple CAA plot
tint=[irf_time([2006 9 27 17 10 0]) irf_time([2006 9 27 17 40 0]) ];
h=irf_plot(4);
irf_plot(h(1),'Data_Velocity_ComponentPerpendicularToMagField__C2_CP__MOMENTS');
irf_plot(h(2),'v_drift_ISR2__C2_CP_EFW_L3_V3D_INERT');
irf_plot(h(3),'E_Vec_xy_ISR2__C2_CP_EFW_L3_E');
irf_plot(h(4),'B_vec_xyz_gse__C2_CP_FGM_5VPS');
irf_zoom(h,'x',tint);
irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
irf_timeaxis(h)
%% Plot panels reading data from other sources than CAA
% copy the following panels into your figure
if 1,   % PANEL: PEACE PADH high resolution pitch C4 (cdf files from QJAS)
  hca=h(i_subplot);i_subplot=i_subplot+1;
  qjas_file='qjas_data_C4_PADH';
  t=irf_cdf_read(qjas_file,'timetags_delta');
  tt=t(:,1);
  dtsampling=0.06;
  dtsampling=0.26;
  psd=cdfread(qjas_file,'VARIABLES','psd');
  theta=cdfread(qjas_file,'VARIABLES','theta');
  phi=cdfread(qjas_file,'VARIABLES','phi');
  level=cdfread(qjas_file,'VARIABLES','level');
  [delmett,ind]=irf_tlim(tt,tint);
  specrec=struct('t',tt(ind),'dt',dtsampling,'p_label','Log PSD [s^6/km^6]');
  psdnew=zeros(length(psd),size(psd{1},1),size(psd{1},2));
  for j=1:length(psd)
    %      flag_measurement=any(psd{j}(:,:),2)*1;
    %      flag_measurement(flag_measurement==0)=NaN;
    %      flag_measurement=repmat(flag_measurement,1,size(psd{j},2));
    psdnew(j,:,:)=psd{j};
    psdnew(psdnew==0)=NaN;
    %      psdnew(j,:,:)=shiftdim(psdnew(j,:,:),1).*flag_measurement;
    %      psdnew(j,:,:)=shiftdim(psdnew(j,:,:),1);
  end
  if 0, % energy spectorgram (integrated over pitch angles)
    specrec.f=log10(res.en);
    specrec.p=res.omni(ind,:);
    specrec.f_label=['Log10 ' res.enlabel];
  elseif 1, % pitch angle spectrogram for given energy
    specrec.f=theta{1};specrec.f_label='Pitch angle';
    enindex=30;level{1}(enindex)
    specrec.f_label=[specrec.f_label '  \newline[E=' num2str(level{1}(enindex),'%6.f') 'eV]'];
    specrec.p=log10(psdnew(ind,:,enindex));
  end
  specrec_C4PADH=specrec;
  irf_spectrogram(hca,specrec);
  %    hold on;
  %    irf_spectrogram(hca,specrec_C43DRH);
  caxis([-2.99 -1.01])
  colormap(jet);
  irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
  set(gca,'ytick',[30 60 90 120 150]);
end
%% Example plots
% see overview of examples under https://sites.google.com/site/andrisvaivads/Andris_Vaivads/cluster/irfu-matlab-examples
open('Example_1.m');
open('Example_2.m');
open('Example_10.m');
%% Cluster data reading from local Uppsala disks
% using c_get_batch
c_get_batch(irf_time([2002 03 04 10 00 00]),30*60,'sp','/home/yuri/caa-data/20020304')
% if time intervals to download are in matrix tint
for j=1:size(tint,1),
  c_get_batch(tint(j,1),tint(j,2)-tint(j,1),'sp',['./' epoch2iso(tint(j,1),1) '-' epoch2iso(tint(j,2),1)]);
end
clear j;

%% Demo of irfu-matlab
% Philosophy
% 1. Maximize time spent on science
% 2. Make data analysis easy and natural
% 3. All instrument data in the same plot
% 4. We want to data analysis ourselves
% 5. Publishing quality figures
% CAA - basic source of data

%____________________________________________________
% Some basic routines and basic usage
irf
irfnotes

% plasma parameters
help irf_plasma_calc
irf_plasma_calc

% cold wave dispersion relation
irf_disp_surf

% plasma wave polarization
irf_plasma_wave_visualization
irf_plasma_wave_visualization('demo','alfven_shear45');

%____________________________________________________
% Example 1
% artificial data
% define time interval
t=irf_time([2008 03 01 10 0 0]):.2:irf_time([2008 03 01 11 0 0]);
t=t(:); % make it column vector
% define one time series
y=exp(0.001*(t-t(1))).*sin(2*pi*t/180);
% define data matrix
data1=[t y];
% plot data
irf_plot(data1)

%____________________________________________________
% Example 2
y=exp(0.001*(t-t(1))).*cos(2*pi*t/180);
data2=[data1 y];
irf_plot(data2)
irf_legend({'X','Y'},[0.02 0.02])

%____________________________________________________
% Example 3
clf;
data3=irf_tappl(data2,'*1.2+2');
h=irf_plot({data2,data3})
ylabel(h(1),'data2');
ylabel(h(2),'data3');
irf_legend(h(1),{'X','Y'},[0.02 0.98],'fontsize',20)

%____________________________________________________
% Example 4
data3=irf_tappl(data2,'*1.2+2');
h=irf_plot({data2,data3},'comp');
ylabel(h(1),'X');
ylabel(h(2),'Y');
irf_legend(h(1),{'data2','data3=data2*1.2+2'},[0.02 0.98],'fontsize',20)
hours=[t (t-t(1))/3600]; % hours from beginning of time interval
irf_timeaxis(h(2),t(1),hours,{'hours'})

%____________________________________________________
% Example 5
clf;
tint1=[irf_time([2008 03 01 10 10 0]) irf_time([2008 03 01 10 20 0])];
tint2=[irf_time([2008 03 01 10 11 0]) irf_time([2008 03 01 10 12 0])];
h=irf_plot(3);
irf_plot(h(1),data2);
ylabel(h(1),'data2');
irf_plot(h(2),data2(:,1:2),'r.','markersize',12);
irf_plot(h(3),data3);
hours=[t (t-t(1))/3600]; % hours from beginning of time interval
irf_zoom(h,'x',tint1);
irf_pl_mark(h(2:3),tint2);
irf_timeaxis(h(end),t(1),hours,{'hours'})
irf_legend(0,'Some additional info',[0,1],'color','r')

%____________________________________________________
% Example 6
% GET REAL DATA
rmdir ~/test;
mkdir ~/test; cd ~/test

tint=[irf_time([2003 02 04 18 0 0]) irf_time([2003 02 04 20 0 0])];
caa_download(tint,'list:*FGM*')
caa_download(tint,'list:*FGM_SPIN*')
caa_download(tint,'*FGM_SPIN*')
% plot
clf
irf_plot('B_vec_xyz_gse__C1_CP_FGM_SPIN');

%____________________________________________________
% Example 7
% Where are satellites?
clf
c_pl_sc_conf_xyz
% to work offline get satellite position data locally
caa_download(tint,'C?_CP_AUX_POSGSE_1M');

%____________________________________________________
% Example 8
% Ions?
caa_download(tint,'list:*CIS*')
caa_download(tint,'list:*CODIF_HS*O*PEF')
caa_download(tint,'*CODIF_O1_1D*PEF')

%plot
clf;
h=irf_plot(2);
irf_plot(h(1),'B_vec_xyz_gse__C4_CP_FGM_SPIN');
irf_plot(h(2),'flux__C4_CP_CIS_CODIF_O1_1D_PEF')

set(h(2),'yscale','log')
irf_zoom(h,'x',tint);
irf_plot_axis_align(2)

%____________________________________________________
% Example 9
% Moments
caa_download(tint,'*CODIF*O*MOMENTS*')

clf;
h=irf_plot(3);
irf_plot(h(1),'B_vec_xyz_gse__C4_CP_FGM_SPIN');
irf_plot(h(2),'flux__C4_CP_CIS_CODIF_O1_1D_PEF')
irf_colormap('space')
irf_plot(h(3),'velocity__C4_CP_CIS_CODIF_HS_O1_MOMENTS')

set(h(2),'yscale','log')
irf_zoom(h,'x',tint);
irf_plot_axis_align(2)

irf_legend(h(3),{'VX','VY','VZ'},[0.02 0.05])

% print
set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% to get bitmap file
print -dpng delme.png
% to get pdf file with no white margins pint to eps and convert after
print -depsc2 -painters delme.eps
% to convert to pdf on the system command line execute some of
% ps2pdf -dEPSFitPage -dEPSCrop delme.eps

%____________________________________________________
% Example 10
% massage of data (MVAR)
[~,~,B1]=c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_SPIN');

clf
irf_plot(B1)

B1=irf_abs(B1);
irf_plot(B1)

% minvar
irf_minvar_gui(B1)

%____________________________________________________
% Example 11
% massage of data (timing)
c_eval('[~,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_SPIN'');');
c_eval('B?=irf_abs(B?);');

c_4_v_gui('B?')

