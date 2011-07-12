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
% more detailed way
% most lines needed to define the size to have best agreement with eps file
set(0,'defaultLineLineWidth', 1.5);
fn=figure(61);clf;
set(fn,'color','white'); % white background for figures (default is grey)
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 24;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xSize sLeft ySize yTop

% set subplots
% specifying position
h(1)=axes('position',[0.65 0.78 0.2 0.2]); % [x y dx dy]
% having all in standard form
n_subplots=8;i_subplot=1;clear h;
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
%% Print the figure as it looks on screen  
set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
print -dpng delme.png
%print -depsc2 delme.eps
%% Add information to figures              
% text and legends
help irf_legend
ht=irf_legend(gca,[mfilename '  ' datestr(now)],[0.02 1.01], 'interpreter','none','fontsize',8);
% labels a),b)...
help irf_pl_number_subplots
% if you want some alternative colorbars, like white in middle, blue negative and red positive
% or if you want to use different colorbars within the same figure 
help irf_colormap
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
%% Cluster data reading from local Uppsala disks
% using c_get_batch
c_get_batch(irf_time([2002 03 04 10 00 00]),30*60,'sp','/home/yuri/caa-data/20020304')
% if time intervals to download are in matrix tint
for j=1:size(tint,1),
    c_get_batch(tint(j,1),tint(j,2)-tint(j,1),'sp',['./' epoch2iso(tint(j,1),1) '-' epoch2iso(tint(j,2),1)]);
end
clear j;
%% Cluster data reading from CAA, data coordinate transformations
% To download from CAA
help caa_download

% To read downloaded CAA data
tint=[irf_time([2006 9 27 17 14 0]) irf_time([2006 9 27 17 26 0])];
if 0, % read s/c position and velocity
    c_eval('[caaR?,~,R?]=c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'');');
    c_eval('[caaV?,~,V?]=c_caa_var_get(''sc_v_xyz_gse__C?_CP_AUX_POSGSE_1M'');');
end
if 0, % read FGM data form all sc
    c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_5VPS'');');
    %    c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');');
    c_eval('B?=irf_abs(B?);');
    c_eval('diB?=c_coord_trans(''GSE'',''ISR2'',B?,''cl_id'',?);');
    c_eval('gsmB?=irf_gse2gsm(B?);');
end
if 0, % read CIS HIA/CODIF velocity moments from available s/c
    c_eval('[caaVCIS?,~,VCIS?]=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',3);
    c_eval('[caaVCISH?,~,VCISH?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',4);
    c_eval('gsmVCIS?=irf_gse2gsm(VCIS?);',3);
    c_eval('gsmVCISH?=irf_gse2gsm(VCISH?);',4);
end
if 0, % read RAPID data
    c_eval('[caaRAPID_J?,~,RAPID_J?]=c_caa_var_get(''Electron_Dif_flux__C?_CP_RAP_ESPCT6'');');
end
if 0, % read EFW data
    c_eval('[caaE?,~,diE?]=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'');');
    c_eval('[caaVps?,~,Vps?]=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'');');
end
%% Cluster data CAA reading and different panel plotting
% example file for creating a figure (run on brain)
cur_dir=pwd;
irf_units;
cd('/share/Cluster/Test/CAA');
tint=[irf_time([2006 9 27 17 10 0]) irf_time([2006 9 27 17 40 0])];
ic=1;
CISinstrument='HIA';
%CISinstrument='CODIF';
if 1, % initialize figure
    fn=figure(61);
    h=irf_plot(8);
    i_subplot=1; % in which subplot is active
    set(fn,'defaultLineLineWidth',1);
end
if 1, % plot figures panels
    if 0, % read FGM data from all sc
        c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');');
        % c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_5VPS'');');
        c_eval('B?=irf_abs(B?);');
        c_eval('gsmB?=irf_gse2gsm(B?);');
    end
    if 1,   % PANEL: C?       FGM B GSM
        hca=irf_panel('C? FGM B GSM');
        c_eval('irf_plot(hca,gsmB?);',ic);
        ylabel(hca,'B [nT] GSM');
        irf_zoom(hca,'y',[-25 25]);
        irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.3])
        irf_legend(hca,{['C' num2str(ic)]},[0.02 0.9],'color','k')
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
    if 1,   % PANEL: C1       CIS HIA/CODIF velocity moment single s/c
        hca=h(i_subplot);i_subplot=i_subplot+1;
        if ic ~=2, % on s/c 2 there is no CIS
            if strcmp(CISinstrument,'HIA')
                dobjname=irf_ssub('C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic);
                varname=irf_ssub('velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS',ic);
            elseif strcmp(CISinstrument,'CODIF')
                dobjname=irf_ssub('C?_CP_CIS_CODIF_HS_H1_MOMENTS',ic);
                varname=irf_ssub('velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS',ic);
            end
            caa_load(dobjname);
            c_eval(['VCIS?=getmat(' dobjname ',''' varname ''');'],ic);
            c_eval('gsmVCIS?=irf_gse2gsm(VCIS?);',ic);
            varunits=eval(['getunits(' dobjname ',''' varname ''')']);
            %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
            disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
            c_eval('irf_plot(hca,gsmVCIS?);',ic);
            ylabel(hca,'V [km/s] GSM');
            set(hca,'ylim',[-199 799]);
            irf_legend(hca,{'V_X','V_Y','V_Z'},[0.02 0.49])
            irf_legend(hca,{['C' num2str(ic)]},[0.02 0.95],'color','k')
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
    if 1,   % PANEL: C1,C3,C4 CIS Vy velocities
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
    if 1,   % PANEL: Pressures, B and CIS HIA/CODIF single s/c
        hca=h(i_subplot);i_subplot=i_subplot+1;
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
    irf_colormap % execute, if necessary with parameter, every time you want to change colormap
    if 1,   % PANEL: CIS HIA/CODIF spectrogram
        hca=h(i_subplot);i_subplot=i_subplot+1;
        if ic~=2,
            if strcmp(CISinstrument,'HIA')
                dobjname=irf_ssub('C?_CP_CIS_HIA_HS_1D_PEF',ic);
                varname=irf_ssub('flux__C?_CP_CIS_HIA_HS_1D_PEF',ic); % HIA
            elseif strcmp(CISinstrument,'CODIF')
                dobjname=irf_ssub('C?_CP_CIS_CODIF_H1_1D_PEF',ic);
                varname=irf_ssub('flux__C?_CP_CIS_CODIF_H1_1D_PEF',ic); % CODIF H+
                %dobjname=irf_ssub('C?_CP_CIS_CODIF_O1_1D_PEF',ic);eval(['caa_load ' dobjname]);varname=irf_ssub('flux__C?_CP_CIS_CODIF_O1_1D_PEF',ic); % CODIF O+
            end
            caa_load(dobjname);
            %varunits=eval(['getunits(' dobjname ',''' varname ''')']);
            varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
            disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
            eval(['plot(hca,' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'');']);
            caxis(hca,[3.9 6.1]);
            irf_colormap;
            set(hca,'yscale','log');
            set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
            ylabel(hca,'E [eV]');
        end
    end
    if 0,   % PANEL: EFW E field in ISR2 reference frame single s/c
        hca=h(i_subplot);i_subplot=i_subplot+1;
        c_eval('irf_plot(hca,diE?)',ic);
        ylabel(hca,'E [mV/m] ISR2');
        irf_zoom(hca,'ylim','smart');
        irf_legend(hca,{'E_X','E_Y'},[0.02 0.49])
        irf_legend(hca,{['C' num2str(ic)]},[0.02 0.95],'color','k')
    end
    if 1,   % PANEL: STAFF spectrogram Bx
        hca=h(i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_STA_PSD',ic);
        caa_load(dobjname);
        varname=irf_ssub('BB_xxyyzz_isr2__C?_CP_STA_PSD',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        varunits='nT^2/Hz';
        disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(hca,' dobjname ',''' varname ''',''ax'',hca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',1);']);
        hold(hca,'on');
        c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic);
        c_eval('flh=irf_plasma_calc(B?,1,0,0,0,''Flh'');',ic);
        irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
        irf_plot(hca,[fce(:,1) fce(:,2)*.5],'-','linewidth',0.2,'color','w');
        irf_plot(hca,[fce(:,1) fce(:,2)*.25],'-','linewidth',0.2,'color','w');
        irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
        caxis(hca,[-10 -7]);
        irf_colormap;
        set(hca,'yscale','lin','ylim',[0 499]);
    end
    if 1,   % PANEL: STAFF spectrogram Ex
        hca=h(i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_STA_PSD',ic);
        caa_load(dobjname);
        varname=irf_ssub('EE_xxyy_isr2__C?_CP_STA_PSD',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        varunits='(mV/m)^2/Hz';
        disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(hca,' dobjname ',''' varname ''',''ax'',hca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',1);']);
        hold(hca,'on');
        c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic);
        c_eval('flh=irf_plasma_calc(B?,1,0,0,0,''Flh'');',ic);
        irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
        irf_plot(hca,[fce(:,1) fce(:,2)*.5],'-','linewidth',0.2,'color','w');
        irf_plot(hca,[fce(:,1) fce(:,2)*.25],'-','linewidth',0.2,'color','w');
        irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
        caxis(hca,[-9 -1]);
        irf_colormap;
        set(hca,'yscale','lin','ylim',[0 599]);
    end
    if 1,   % PANEL: RAPID spectrogram
        hca=h(i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_RAP_ESPCT6',ic);
        caa_load(dobjname);
        varname=irf_ssub('Electron_Dif_flux__C?_CP_RAP_ESPCT6',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
        disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(hca,' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'');']);
        caxis(hca,[0.51 4.49]);
        irf_colormap;
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
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
    if 1,   % PANEL: PEACE PITCH_SPIN_DEFlux spectrogram omni
        hca=h(i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        caa_load(dobjname);
        varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
        disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(hca,' dobjname ',''' varname ''',''sum_dim1'',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'');']);
        caxis(hca,[5.8 7.6]);
        irf_colormap;
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E [eV]');
    end
    if 0,   % PANEL: PEACE PITCH_SPIN_DEFlux spectrogram angles
        hca=h(i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        caa_load(dobjname);
        varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
        disp(['PANEL: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(hca,' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',15);']);
        caxis(hca,[5.8 7.6]);
        irf_colormap;
        set(hca,'yscale','lin');
        set(hca,'ytick',[0 45 90 135 180]);
        ylabel(hca,'\Theta [deg]');
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
    % subspin resolution panels
    if 0,   % PANEL: PEACE PEA_PITCH_3DRH_PSD high res
        hca=h(i_subplot);i_subplot=i_subplot+1;
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
    irf_plot_axis_align
    irf_zoom(h,'x',tint);
    irf_legend(h(1),'Figure reference',[0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);    irf_plot_axis_align
    irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
    irf_timeaxis(h);
end
cd(cur_dir)
%print -depsc2 -painters fig/vaivads2011a_fig2.eps
%c_eval('print -dpng -painters fig/vaivads2011a_Figure_2z_C?.png',ic);
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


%%  Figure  3
% top - spatial scale
% 1. 4 sc B
% 2. density, CIS + EFW
% 3. E Ex,Ey
% 4. RAPID pitch angle
% 5. Bx waves


%% Figure s/c configuration 
