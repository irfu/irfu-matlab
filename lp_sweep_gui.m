function lp_sweep_gui(action)
%LP_SWEEP_GUI interactively work with sweeps
%
% LP_SWEEP
%
% You can access the results from the figures data 'userdata'.
% data=get(gcf,'userdata');ud=data.ud
%
% ud.I - current
% ud.U - voltage
% ud.dUdI - resistance
% ud.... - many other parameters
%
irf_units;
persistent message;
if isempty(message), % run only the first time during the session
    message='You can anytime access all the results from get(gcf,''userdata'').';
    disp(message);
end
if      nargin == 0, action='initialize';end
switch action,
    case 'initialize'
        %% default values
        ud=struct();
        ud.U_string='-10:.1:50';
        ud.R_sun=1; % distance in AU
        ud.UV_factor=1;
        % probe
        ud.probe.type='spherical';
        ud.probe.surface='themis'; % should be one of options in lp_photocurrent
        ud.probe.radius=4; % in cm
        ud.probe.length=4; % in cm
        ud.probe.total_vs_sunlit_area=4;
        ud.probe.total_area=(ud.probe.radius*.01)^2*4*pi;
        % s/c
        ud.sc_radius=1; % [m]
        ud.flag_use_sc=0; % 0-not use, 1- use sc
        ud.sc_probe_refpot_as_fraction_of_scpot=0.25; % reference potential at probe
        ud.sc.number_of_probes=4;
        ud.sc.illuminated_area=1; % illuminated cross section generating photoelectrons
        ud.sc.total_area=4;       % total s/c area collecting photoelectrons
        % plasma
        ud.m_amu1=1;
        ud.m_amu2=16;
        ud.m2=0; % relative fraction of 2nd species
        ud.n=1; % [cc]
        ud.Ti=1;
        ud.Te=1;
        ud.V_SC=0;
        
        %% initialize figure
        set(0,'defaultLineLineWidth', 1.5);
        fn=figure(63);
        clf reset;
        clear h;
        set(fn,'color','white'); % white background for figures (default is grey)
        set(gcf,'PaperUnits','centimeters')
        set(gcf,'defaultAxesFontSize',14);
        set(gcf,'defaultTextFontSize',14);
        set(gcf,'defaultAxesFontUnits','pixels');
        set(gcf,'defaultTextFontUnits','pixels');
        xSize = 13; ySize = 13;
        xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
        set(gcf,'Position',[100 300 xSize*50 ySize*50])
        set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
        clear xSize sLeft ySize yTop
        %        set(fn,    'windowbuttondownfcn', 'irf_minvar_gui(''ax'')');zoom off;
        ud.h(1)=axes('position',[0.1 0.25 0.5 0.3]); % [x y dx dy]
        ud.h(2)=axes('position',[0.1 0.65 0.5 0.3]); % [x y dx dy]
        ud.h(3)=axes('position',[0.1 0.0 0.5 0.13]); % [x y dx dy]
        linkaxes(ud.h(1:2),'x');
        axis(ud.h(3),'off');
        ud.ht=text(0,1,'','parent',ud.h(3));
        set(fn,'userdata',ud);
        
        %% initialize probe menu
        hp = uipanel('Title','Probe','FontSize',12,'BackgroundColor',[1 0.95 1],'Position',[.7 .0 .3 .42]);
        inp.update                           = uicontrol('Parent',hp,'String','Update',                'Position',[0   0 60 30],'Callback','lp_sweep_gui(''update'')');
        inp.reset                            = uicontrol('Parent',hp,'String','Reset',                 'Position',[70  0 60 30],'callback', 'lp_sweep_gui(''initialize'')');
        inp.U_text                           = uicontrol('Parent',hp,'String','U [V]',                 'Position',[0   25 60 30]);
        inp.U_value                          = uicontrol('Parent',hp,'String',ud.U_string,             'Position',[70  25 100 30],'style','edit','backgroundcolor','white');
        inp.Rsun_text                        = uicontrol('Parent',hp,'String','Rsun [AU]',             'Position',[0   50 60 30]);
        inp.Rsun_value                       = uicontrol('Parent',hp,'String',num2str(ud.R_sun),       'Position',[70  50 100 30],'style','edit','backgroundcolor','white');
        inp.UV_factor_text                   = uicontrol('Parent',hp,'String','UV factor',             'Position',[0   75 60 30]);
        inp.UV_factor_value                  = uicontrol('Parent',hp,'String',num2str(ud.UV_factor),   'Position',[70  75 100 30],'style','edit','backgroundcolor','white');
        inp.probe_radius_text                = uicontrol('Parent',hp,'String','probe radius [cm]',     'Position',[0   125 90 30]);
        inp.probe_radius_value               = uicontrol('Parent',hp,'String',num2str(ud.probe.radius),'Position',[100 125 70 30],'style','edit','backgroundcolor','white');
        inp.probe_length_text                = uicontrol('Parent',hp,'String','probe length [cm]',     'Position',[0   150 90 30]);
        inp.probe_length_value               = uicontrol('Parent',hp,'String','','style','text',       'Position',[100 150 70 30],'backgroundcolor','white');
        inp.probe_total_vs_sunlit_area_text  = uicontrol('Parent',hp,'String','total/sunlit area',     'Position',[0   175 90 30]);
        inp.probe_total_vs_sunlit_area_value = uicontrol('Parent',hp,'String',num2str(ud.probe.total_vs_sunlit_area),'style','text','Position',[100 175 70 30],'backgroundcolor','white');
        inp.probe_type                       = uicontrol('Parent',hp,'String','spherical probe|cylindrical probe|specify probe area','style','popup','Position',[2 225 150 30],'backgroundcolor','white','Callback', @setprobetype);
        inp.probe_surface                    = uicontrol('Parent',hp,'String','probe_surface|themis|cassini',         'Position',[2 200 150 30],'style','popup','backgroundcolor','white','Callback', @setprobesurface);
        
        %% initialize s/c menu
        hsc = uipanel('Title','Spacecraft','FontSize',12,'BackgroundColor',[.95 1 1],'Position',[.7 .43 .3 .27]);
        inp.flag_sc                                    = uicontrol('Parent',hsc,'style','radio','String','Model spacecraft','Value',0,             'Position',[0 125 120 25]);
        inp.sc_example                                 = uicontrol('Parent',hsc,'String','Example spacecraft|Cluster|Solar Orbiter|THEMIS|Cassini','Position',[0 100 150 25],'style','popup','backgroundcolor','white','Callback', @setscexample);
        inp.sc.illuminated_area_text                   = uicontrol('Parent',hsc,'String','Illuminated area [m2]',                                  'Position',[0 50 120 25]);
        inp.sc.illuminated_area_value                  = uicontrol('Parent',hsc,'String',num2str(ud.sc.illuminated_area),'style','edit',           'Position',[120 50 50 25],'backgroundcolor','white');
        inp.sc.total_area_text                         = uicontrol('Parent',hsc,'String','Total area [m2]',                                        'Position',[0 75 120 25]);
        inp.sc.total_area_value                        = uicontrol('Parent',hsc,'String',num2str(ud.sc.total_area),'style','edit',                 'Position',[120 75 50 25],'backgroundcolor','white');
        inp.sc_probe_refpot_as_fraction_of_scpot_text  = uicontrol('Parent',hsc,'String','Probe refpot/scpot',                                     'Position',[0 25 120 25]);
        inp.sc_probe_refpot_as_fraction_of_scpot_value = uicontrol('Parent',hsc,'String',num2str(ud.sc_probe_refpot_as_fraction_of_scpot),         'Position',[120 25 50 25],'style','edit','backgroundcolor','white');
        inp.sc_number_of_probes_text                   = uicontrol('Parent',hsc,'String','Number of probes',                                       'Position',[0 0 120 25]);
        inp.sc_number_of_probes_value                  = uicontrol('Parent',hsc,'String',num2str(ud.sc.number_of_probes),                          'Position',[120 0 50 25],'style','edit','backgroundcolor','white');
        
        %% initialize plasma menu
        hpl= uipanel('Title','Plasma','FontSize',12,'BackgroundColor',[1 1 .95],'Position',[.7 .7 .3 .25]);
        inp.n         = uicontrol('Parent',hpl,'String','Ne [cc]','Position',[0 0 80 25]);
        inp.n_value   = uicontrol('Parent',hpl,'String',num2str(ud.n),'style','edit','Position',[80 0 90 25],'backgroundcolor','white');
        inp.T         = uicontrol('Parent',hpl,'String','T [eV]','Position',[0 25 80 25]);
        inp.T_value   = uicontrol('Parent',hpl,'String','1 1','style','edit','Position',[80 25 90 25],'backgroundcolor','white');
        inp.m         = uicontrol('Parent',hpl,'String','m [mp],0=me','Position',[0 50 80 25]);
        inp.m_value   = uicontrol('Parent',hpl,'String','0 1','style','edit','Position',[80 50 90 25],'backgroundcolor','white');
        inp.q         = uicontrol('Parent',hpl,'String','q [e]','Position',[0 75 80 25]);
        inp.q_value   = uicontrol('Parent',hpl,'String','-1 1','style','edit','Position',[80 75 90 25],'backgroundcolor','white');
        inp.vsc       = uicontrol('Parent',hpl,'String','Vsc [m/s]','Position',[0 100 80 25]);
        inp.vsc_value = uicontrol('Parent',hpl,'String','0','style','edit','Position',[80 100 90 25],'backgroundcolor','white');
        
        ud.inp=inp;
        set(gcf,'userdata',ud);
        %
        lp_sweep_gui('update');
    case 'update'
        disp('Updating');
        ud=get(gcf,'userdata');
        %% get input parameters
        inp=ud.inp;
        ud.U=eval(get(inp.U_value,'string'));
        ud.UV_factor=str2double(get(inp.UV_factor_value,'string'));
        ud.R_sun=str2double(get(inp.Rsun_value,'string'));
        ud.probe.radius=str2double(get(inp.probe_radius_value,'string'));
        ud.probe.length=str2double(get(inp.probe_length_value,'string'));
        ud.n  =eval(['[' get(inp.n_value,'string')   ']' ]);
        ud.T  =eval(['[' get(inp.T_value,'string')   ']' ]);
        ud.m  =eval(['[' get(inp.m_value,'string')   ']' ]);
        ud.q  =eval(['[' get(inp.q_value,'string')   ']' ]);
        ud.vsc=eval(['[' get(inp.vsc_value,'string') ']' ]);       
        % if scflag read in sc parameters
        ud.flag_use_sc=get(inp.flag_sc,'Value');
        if ud.flag_use_sc,
            ud.probe_refpot_as_fraction_of_scpot=str2double(get(inp.sc_probe_refpot_as_fraction_of_scpot_value,'string'));
            ud.sc.number_of_probes=str2double(get(inp.sc_number_of_probes_value,'string')); % in cm
            ud.sc.cross_section_area=str2double(get(inp.sc.illuminated_area_value,'string'));
            ud.sc.total_area=str2double(get(inp.sc.total_area_value,'string'));
            ud.sc.type='spherical';
            ud.sc.surface='default';
        end
        %% calculate IU curves
        Upot=ud.U(:);
        probe=ud.probe;
        switch ud.probe.type
            case 'spherical'
                probe.cross_section_area=pi*(ud.probe.radius*.01)^2;
                probe.total_area=4*probe.cross_section_area;
                probe.capacitance=4*pi*Units.eps0*ud.probe.radius*.01;
            case 'cylindrical'
                probe.cross_section_area=2*ud.probe.radius*ud.probe.length*0.01^2;
                probe.total_area=pi*probe.cross_section_area;
                probe.length=str2double(get(inp.probe_length_value,'string'));
                probe.capacitance=2*pi*Units.eps0*ud.probe.length*0.01/log(ud.probe.length/ud.probe.radius); % assuming length >> radius
            case 'arbitrary'
                probe.type='spherical';
        end
        ud.probe=probe;
        [J_probe,J_photo,J_plasma]=lp_probe_current(probe,Upot,ud.R_sun,ud.UV_factor,ud);
        dUdI=gradient(Upot,J_probe);
        ud.I=J_probe;
        ud.dUdI=dUdI;
        % if scflag then calculate s/c IU curve
        if ud.flag_use_sc,
            J_sc=lp_probe_current(ud.sc,Upot,ud.R_sun,ud.UV_factor,ud);
            ud.I_sc=J_sc;
            ud.dUdI_sc=gradient(Upot,J_sc);
        end
        % if scflag calculate probe to sc IU curve
        if ud.flag_use_sc,
            Iprobe=min(J_probe):.01*(max(J_probe)-min(J_probe)):max(J_probe);
            Iprobe=Iprobe(:);
            Isat=-ud.sc.number_of_probes*Iprobe;
            Usatsweep=interp1(J_sc,Upot,Isat); % floating potential of sc during sweep
            % plasma current with UV factor zero
            J_probe_plasma=lp_probe_current(probe,Upot,ud.R_sun,0.00000000,ud);
            Uproberefsweep=ud.probe_refpot_as_fraction_of_scpot*Usatsweep; % reference potential around probe
            Uprobe2plasma=zeros(size(Iprobe)); % initialize
            for ii=1:numel(Iprobe),
                % photoelectron current with plasma current zero
                [~,J_probe_photo]=lp_probe_current(probe,Upot-Uproberefsweep(ii),ud.R_sun,ud.UV_factor,ud);
                J_probe=J_probe_plasma+J_probe_photo;
                Uprobe2plasma(ii)=interp1(J_probe,Upot,Iprobe(ii));
            end
            
            Uprobe2sc        =Uprobe2plasma-Usatsweep;
            Uprobe2refpot    =Uprobe2plasma-Uproberefsweep;
            dUdI_probe2plasma=gradient(Uprobe2plasma,Iprobe);
            dUdI_probe2sc    =gradient(Uprobe2sc,Iprobe);
            dUdI_probe2refpot=gradient(Uprobe2refpot,Iprobe);
            dUdI             =dUdI_probe2refpot;
            Upot             =Uprobe2refpot;
            J_probe          =Iprobe;
        end
        %% plot IU curve
        info_txt='';
        h=ud.h;
        plot(h(1),Upot,J_probe*1e6,'k');
        grid(h(1),'on');
        xlabel(h(1),'U [V]');
        ylabel(h(1),'I [\mu A]');
        if ud.flag_use_sc, % add probe potential wrt s/c and plasma
            hold(h(1),'on');
            plot(h(1),Uprobe2plasma,Iprobe*1e6,'r','linewidth',1.5);
            plot(h(1),Uprobe2sc,Iprobe*1e6,'b','linewidth',1.5);
            plot(h(1),Usatsweep,Iprobe*1e6,'color',[0.5 0 0.5],'linewidth',1);
            hold(h(1),'off');
            irf_legend(h(1),'probe to plasma',[0.02 0.98],'color','r');
            irf_legend(h(1),'probe to reference',[0.02 0.88],'color','k');
            irf_legend(h(1),'probe to s/c',[0.02 0.78],'color','b');
            irf_legend(h(1),'Satellite potential',[0.98 0.03],'color',[0.5 0 0.5]);
        else % add photoelectron and photoelectron currents
            hold(h(1),'on');
            plot(h(1),Upot,J_photo*1e6,'r','linewidth',0.5);
            irf_legend(h(1),'    total',      [0.98 0.03],'color','k');
            irf_legend(h(1),' photoelectrons',[0.98 0.13],'color','r');
            clr=[0.5 0 0; 0 0.5 0; 0 0 0.5];
            for ii=1:length(J_plasma),
                plot(h(1),Upot,J_plasma{ii}*1e6,'linewidth',.5,'color',clr(:,ii));
                irf_legend(h(1),['plasma ' num2str(ii)],[0.98 0.13+ii*0.1],'color',clr(:,ii));
            end
            hold(h(1),'off');            
        end
        
        plot(h(2),Upot,dUdI,'k');
        grid(h(2),'on');xlabel(h(2),'U [V]');
        ylabel(h(2),'dU/dI [\Omega]');
        if ud.flag_use_sc, % add probe resistance wrt plasma and s/c
            hold(h(2),'on');
            plot(h(2),Uprobe2plasma,dUdI_probe2plasma,'r','linewidth',1.5);
            plot(h(2),Uprobe2sc,dUdI_probe2sc,'b','linewidth',1.5);
            hold(h(2),'off');
        end
        axis(h(2),'tight');
        Rmin = min(abs(dUdI)); % minimum resistance
        fcr=1/2/pi/Rmin/probe.capacitance;
        disp(['Rmin=' num2str(Rmin,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, f_{CR}=' num2str(fcr,3) 'Hz.']);
        if ud.flag_use_sc,
            info_txt=[info_txt '\newline probe to plasma Rmin=' num2str(min(abs(dUdI_probe2plasma)),3) ' Ohm'];
            info_txt=[info_txt '\newline probe to reference Rmin=' num2str(min(abs(dUdI)),3) ' Ohm'];
            info_txt=[info_txt '\newline probe to spacecraft Rmin=' num2str(min(abs(dUdI_probe2sc)),3) ' Ohm'];
        else
            info_txt=[info_txt '\newline Rmin=' num2str(Rmin,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, f_{CR}=' num2str(fcr,3) 'Hz.'];
        end
        if min(J_probe)<0 && max(J_probe)>0, % display information on Ufloat
            Ufloat=interp1(J_probe,Upot,0); % floating potential
            ii=isfinite(Upot);Rfloat=interp1(Upot(ii),dUdI(ii),Ufloat);
            info_txt=[info_txt '\newline Ufloat=' num2str(Ufloat,3) 'V, R at Ufloat R= ' num2str(Rfloat,3) ' Ohm'];
            disp(['Ufloat=' num2str(Ufloat,3) ' V, R at Ufloat R=' num2str(Rfloat,3) ' Ohm']);
        end
        if ud.UV_factor>0,  % display photoelectron saturation current
            info_txt=[info_txt '\newline photo e- Io = ' num2str(ud.UV_factor*lp_photocurrent(1,-1,ud.R_sun,ud.probe.surface)*1e6,3) '[\mu A/m^2]'];
            if ud.R_sun~=1,
                info_txt=[info_txt '  (' num2str(ud.UV_factor*lp_photocurrent(1,-1,1,ud.probe.surface)*1e6,3) ' \mu A/m^2 at 1 AU)'];
            end
        end
        set(h(2),'yscale','log')
        
        
        axis(h(3),'off');
        set(ud.ht,'string',info_txt);
        set(gcf,'userdata',ud);
end
end

function setprobetype(hObj,event) %#ok<INUSD>
% Called when user activates popup menu of minvar method
val = get(hObj,'Value');
data=get(gcf,'userdata');
if val ==1
    data.probe.type='spherical';
    set(data.inp.probe_length_value,'style','text','string','');
elseif val == 2
    data.probe_type='cylindrical';
    set(data.inp.probe_length_value,'style','edit','string',num2str(data.probe.length));
elseif val == 3
    data.probe.type='spherical';
end
set(gcf,'userdata',data);
end
function setprobesurface(hObj,event) %#ok<INUSD>
% Called when user activates popup menu of minvar method
val = get(hObj,'Value');
data=get(gcf,'userdata');
if val ==1
    data.probe.surface='themis';
elseif val == 2
    data.probe.surface='themis';
elseif val == 3
    data.probe.surface='cassini';
end
set(gcf,'userdata',data);
end
function setscexample(hObj,event) %#ok<INUSD>
% Called when user activates popup menu of minvar method
val = get(hObj,'Value');
data=get(gcf,'userdata');
if val ==1 % do nothing, shows in menu 'Example spacecraft'
elseif val ==2, % Cluster s/c diameter 2.9, height 1.3m
    data.probe.type='spherical';
    set(data.inp.probe_type,'Value',1);
    data.probe.surface='themis';
    set(data.inp.probe_surface,'Value',2);
    set(data.inp.probe_length_value,'style','text','string','');
    set(data.inp.probe_radius_value,'string','4');
    set(data.inp.sc.illuminated_area_value,'string','3.87');
    set(data.inp.sc.total_area_value,'string','25.66');
    set(data.inp.sc_probe_refpot_as_fraction_of_scpot_value,'string','.2');
    set(data.inp.sc_number_of_probes_value,'string','4');
    set(data.inp.Rsun_value,'string','1');
    set(data.inp.probe_total_vs_sunlit_area_value,'string','4');
    data.probe_total_vs_sunlit_area=4;
    set(data.inp.n_value,'string','1');
    set(data.inp.T_value,'string','100 500');
elseif val == 3 % Solar Orbiter with solar panel back side
    data.probe.type='cylindrical';
    set(data.inp.probe_type,'Value',2);
    data.probe.surface='themis';
    set(data.inp.probe_surface,'Value',2);
    set(data.inp.probe_radius_value,'string','0.575');
    set(data.inp.probe_length_value,'style','edit','string','500');
    set(data.inp.sc.illuminated_area_value,'string','5.95');
    set(data.inp.sc.total_area_value,'string','28.11');
    set(data.inp.probe_total_vs_sunlit_area_value,'string',num2str(pi,4));
    set(data.inp.sc_probe_refpot_as_fraction_of_scpot_value,'string','.2');
    set(data.inp.sc_number_of_probes_value,'string','3');
    data.probe_total_vs_sunlit_area=pi;
    set(data.inp.Rsun_value,'string','0.28');
    set(data.inp.n_value,'string','100');
    set(data.inp.T_value,'string','20 40');
elseif val == 4 % THEMIS
    data.probe.type='spherical';
    set(data.inp.probe_type,'Value',1);
    data.probe.surface='themis';
    set(data.inp.probe_surface,'Value',2);
    set(data.inp.probe_length_value,'style','text','string','');
    set(data.inp.probe_radius_value,'string','4');
    set(data.inp.sc.illuminated_area_value,'string','0.71');
    set(data.inp.sc.total_area_value,'string','4.2');
    set(data.inp.sc_probe_refpot_as_fraction_of_scpot_value,'string','.2');
    set(data.inp.sc_number_of_probes_value,'string','4');
    set(data.inp.Rsun_value,'string','1');
    set(data.inp.probe_total_vs_sunlit_area_value,'string','4');
    data.probe_total_vs_sunlit_area=4;
    set(data.inp.n_value,'string','1');
    set(data.inp.T_value,'string','100 500');
elseif val == 5 % Cassini, sensor - 50mm sphere on 10.9 cm stub (diameter 6.35 mm) (efficient radius 56.5 taking into account stub)
    % sc - >6.7 metres high and >4 metres wide
    data.probe.type='spherical';
    set(data.inp.probe_type,'Value',1);
    data.probe.surface='cassini';
    set(data.inp.probe_surface,'Value',3);
    set(data.inp.probe_length_value,'style','text','string','');
    set(data.inp.probe_radius_value,'string','2.5');
    set(data.inp.sc.illuminated_area_value,'string','27');
    set(data.inp.sc.total_area_value,'string','140');
    set(data.inp.sc_probe_refpot_as_fraction_of_scpot_value,'string','.2');
    set(data.inp.sc_number_of_probes_value,'string','1');
    set(data.inp.Rsun_value,'string','9.5');
    set(data.inp.probe_total_vs_sunlit_area_value,'string','4');
    data.probe.total_vs_sunlit_area=4;
    set(data.inp.n_value,'string','0.1');
    set(data.inp.T_value,'string','1 1');
end
set(gcf,'userdata',data);
end
