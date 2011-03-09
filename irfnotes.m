% IRFNOTES File with different common examples how to use irf routines
% enable code folding (including cells) to fast find your necessary examples
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
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
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
% example file for creating figure (run on brain)
cur_dir=pwd;
irf_units;
cd('/share/Cluster/Test/CAA');
tint=[toepoch([2006 9 27 17 10 0]) toepoch([2006 9 27 17 40 0])];
ic=1;
CISinstrument='HIA';
%CISinstrument='CODIF';
if 1,% initialize figure
    set(0,'defaultLineLineWidth', 1.5);
    fn=figure(61); clf;clear h;
    set(fn,'color','white'); % white background for figures (default is grey)
    set(gcf,'PaperUnits','centimeters')
    xSize = 12; ySize = 24;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 10 xSize*50 ySize*50])
    % set subplots
    % specifying position
    %h(1)=axes('position',[0.65 0.78 0.2 0.2]); % [x y dx dy]
    % having all in standard form
    clear xSize sLeft ySize yTop
end
n_subplots=8;i_subplot=1;
if 1, % plot figures panels
    if 1,   % PANEL: FGM B single s/c
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_FGM_5VPS',ic);
        caa_load(dobjname);
        varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_5VPS',ic);
        c_eval(['B?=getmat(' dobjname ',''' varname ''');'],ic);
        c_eval('B?=irf_abs(B?);',ic);
        c_eval('gsmB?=irf_gse2gsm(B?);',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        c_eval('irf_plot gsmB?',ic);
        ylabel('B [nT] GSM');set(gca,'ylim',[-24.9 24.9]);
        irf_legend(gca,{'B_X','B_Y','B_Z','B'},[0.02 0.3])
        irf_legend(gca,{['C' num2str(ic)]},[0.02 0.9],'color','k')
    end
    irf_legend(gca,'Figure reference',[0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);    irf_plot_axis_align
    if 1,   % PANEL: CIS HIA/CODIF velocity moment single s/c
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
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
            disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
            c_eval('irf_plot gsmVCIS?',ic);
            ylabel('V [km/s] GSM');set(gca,'ylim',[-199 799]);
            irf_legend(gca,{'V_X','V_Y','V_Z'},[0.02 0.49])
            irf_legend(gca,{['C' num2str(ic)]},[0.02 0.95],'color','k')
        end
    end
    if 1,   % PANEL: Pressures, B and CIS HIA/CODIF single s/c
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
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
            c_eval(['VelocityCIS?=irf_abs(VelocityCIS?);'],ic); %
            c_eval(['DynamicPressureCIS?=irf_multiply(0.5*1e6*Units.mp*1e6*1e9,DensityCIS?,1,VelocityCIS?(:,[1 5]),2);'],ic); % in nPa, assumes H+
            c_eval(['PressureB?=irf_multiply(1e-18*0.5/Units.mu0*1e9,B?(:,[1 5]),2,B?(:,[1 5]),0);'],ic); % Pressure in nPa
            c_eval(['PressureTotal?=irf_add(1,PressureCIS?,1,PressureB?);'],ic); % Pressure in nPa
            varunits=eval(['getunits(' dobjname ',''' varname ''')']);
            %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
            disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
            c_eval('irf_plot({PressureB?,PressureCIS?,DynamicPressureCIS?,PressureTotal?},''comp'');',ic);
            ylabel('P [nPa]');set(gca,'ylim',[0 0.29]);
            irf_legend(gca,{'P_B','P_i','P_{kin}'},[0.02 0.3])
            irf_legend(gca,{['C' num2str(ic)]},[0.02 0.9],'color','k')
        end
    end
    if 1,   % PANEL: CIS HIA/CODIF spectrogram
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
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
            disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
            eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'');']);
            caxis([3.9 6.1]);
            irf_colormap;
            set(gca,'yscale','log'); set(gca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
            ylabel('E [eV]');
        end
    end
    if 0,   % PANEL: EFW E field in ISR2 reference frame single s/c
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_EFW_L2_E',ic);
        varname=irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L2_E',ic);
        caa_load(dobjname);
        c_eval(['diE?=getmat(' dobjname ',''' varname ''');'],ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
        disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        c_eval('irf_plot diE?',ic);
        ylabel('E [mV/m] ISR2');set(gca,'ylim',[-10 10]);
        irf_legend(gca,{'E_X','E_Y'},[0.02 0.49])
        irf_legend(gca,{['C' num2str(ic)]},[0.02 0.95],'color','k')
    end
    if 1,   % PANEL: RAPID spectrogram
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_RAP_ESPCT6',ic);
        caa_load(dobjname);
        varname=irf_ssub('Electron_Dif_flux__C?_CP_RAP_ESPCT6',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
        disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'');']);
        caxis([0.51 4.49]);
        irf_colormap;
        set(gca,'yscale','log'); set(gca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
    end
    if 0,   % PANEL: RAPID spectrogram parallel
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
        caa_load(dobjname);
        varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
        l3dd_caa=eval(['getv(' dobjname ',''' varname ''')']);
        l3dd=eval(['getmat(' dobjname ',''' varname ''')']);
        l3dd_energies=eval(['getv(' dobjname ',''' l3dd_caa.DEPEND_1 ''')']);
        disp(['SUBPLOT: C' num2str(ic) ' RAPID anisotropy']);
        disp(['dobj:' dobjname ]);disp([' var:' varname]);
        eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',1);']);
        caxis([1.1 4.3]);
        irf_colormap;
        set(gca,'yscale','log'); set(gca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5]);
    end
    if 0,   % PANEL: RAPID spectrogram perpendicular
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
        caa_load(dobjname);
        varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
        l3dd_caa=eval(['getv(' dobjname ',''' varname ''')']);
        l3dd=eval(['getmat(' dobjname ',''' varname ''')']);
        l3dd_energies=eval(['getv(' dobjname ',''' l3dd_caa.DEPEND_1 ''')']);
        disp(['SUBPLOT: C' num2str(ic) ' RAPID anisotropy']);
        disp(['dobj:' dobjname ]);disp([' var:' varname]);
        eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',5);']);
        caxis([1.1 4.3]);
        irf_colormap;
        set(gca,'yscale','log'); set(gca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5]);
    end
    if 0,   % PANEL: RAPID spectrogram pitch angle
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
        caa_load(dobjname);
        varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
        disp(['SUBPLOT: C' num2str(ic) ' RAPID anisotropy']);
        disp(['dobj:' dobjname ]);disp([' var:' varname]);
        eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp_dim1'',''comp'',1);']);
        caxis([1.1 4.3]);
        irf_colormap;
        set(gca,'yscale','lin'); 
        set(gca,'ytick',[0 45 90 135 180]);
        ylabel(gca,'Pitch ang. [deg]');
    end
    if 0,   % PANEL: RAPID spectrogram anisotropy
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_RAP_PAD_L3DD',ic);
        caa_load(dobjname);
        varname=irf_ssub('PAD_Electron_L_Dif_flux__C?_CP_RAP_PAD_L3DD',ic);
        l3dd_caa=eval(['getv(' dobjname ',''' varname ''')']);
        l3dd=eval(['getmat(' dobjname ',''' varname ''')']);
        l3dd_energies=eval(['getv(' dobjname ',''' l3dd_caa.DEPEND_1 ''')']);
        disp(['SUBPLOT: C' num2str(ic) ' RAPID anisotropy']);
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
        caa_spectrogram(gca,specrec);
        colorbar;
        caxis([-2 2]);
        irf_colormap('poynting');
        set(gca,'yscale','log'); set(gca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
    end
    if 1,   % PANEL: PEACE spectrogram omni
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        caa_load(dobjname);
        varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
        disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(' dobjname ',''' varname ''',''sum_dim1'',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'');']);
        caxis([5.8 7.6]); 
        irf_colormap;
        set(gca,'yscale','log'); set(gca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel('E [eV]');
    end
    if 0,   % PANEL: PEACE spectrogram angles
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        caa_load(dobjname);
        varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
        disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',15);']);
        caxis([5.8 7.6]); 
        irf_colormap;
        set(gca,'yscale','lin'); 
        set(gca,'ytick',[0 45 90 135 180]);
        ylabel(gca,'\Theta [deg]');
    end
    if 1,   % STAFF: spectrogram Bx
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_STA_PSD',ic);
        caa_load(dobjname);
        varname=irf_ssub('BB_xxyyzz_isr2__C?_CP_STA_PSD',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        varunits='nT^2/Hz';
        disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',1);']);
        hold on;
        c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic);
        irf_plot(fce,'-','linewidth',0.2,'color','k');
        caxis([-10 -7]); 
        irf_colormap;
        set(gca,'yscale','lin','ylim',[0 499]); 
    end
    if 1,   % STAFF: spectrogram Ex
        h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
        dobjname=irf_ssub('C?_CP_STA_PSD',ic);
        caa_load(dobjname);
        varname=irf_ssub('EE_xxyy_isr2__C?_CP_STA_PSD',ic);
        varunits=eval(['getunits(' dobjname ',''' varname ''')']);
        varunits='(mV/m)^2/Hz';
        disp(['SUBPLOT: C' num2str(ic)]);disp(['dobj:' dobjname ]);disp([' var:' varname]);disp(['varunits: ' varunits]);
        eval(['plot(' dobjname ',''' varname ''',''ax'',gca,''colorbarlabel'',''' varunits ''',''fitcolorbarlabel'',''comp'',1);']);
        hold on;
        c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic);
        irf_plot(fce,'-','linewidth',0.2,'color','k');
        caxis([-9 -1]); 
        irf_colormap;
        set(gca,'yscale','lin','ylim',[0 599]); 
    end
    irf_plot_axis_align
    irf_zoom(tint,'x',h);
    add_timeaxis(h);
end
cd(cur_dir)
%print -depsc2 -painters fig/vaivads2011a_fig2.eps
%c_eval('print -dpng -painters fig/vaivads2011a_Figure_2z_C?.png',ic);
%% Colors/colorbars in matlab
% colorbar with white in middle, blue negative and red positive
help irf_colromap
%% Other
