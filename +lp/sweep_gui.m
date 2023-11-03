function sweep_gui(action)
%LP.SWEEP_GUI interactively work with sweeps
%
% LP.SWEEP_GUI
%
% You can access the results from the figures data 'userdata'.
% data=get(gcf,'userdata');ud=data.ud
%
% ud.I - current
% ud.U - voltage
% ud.dUdI - resistance
% ud.... - many other parameters
%
Units=irf_units;
persistent message;
if isempty(message) % run only the first time during the session
  message='You can anytime access all the results from get(gcf,''userdata'').';
  disp(message);
end
if      nargin == 0, action='initialize';end
switch action
  case 'initialize'
    %% default values
    ud=struct();
    ud.U_string='-10:.1:50';
    ud.R_sun=1; % distance in AU
    ud.UV_factor=1;
    % probe
    ud.probe.type='spherical';
    ud.probe.surface='themis'; % should be one of options in lp.photocurrent
    ud.probe.radius=4; % in cm
    ud.probe.length=4; % in cm
    ud.probe.total_vs_sunlit_area=4;
    ud.probe.total_area=(ud.probe.radius*.01)^2*4*pi;
    % s/c
    ud.sc_radius=1; % [m]
    ud.flag_use_sc=0; % 0-not use, 1- use sc
    ud.sc.probe_refpot_as_fraction_of_scpot=0.25; % reference potential at probe
    ud.sc.number_of_probes=4;
    ud.sc.probe_distance_to_spacecraft=44;
    ud.sc.sunlit_area=1;      % sunlit cross section generating photoelectrons
    ud.sc.total_area=4;       % total s/c area collecting photoelectrons
    % plasma default values are specified in plasma menu (can be moved here
    % if deemed necessary
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
    xSize = 13; ySize = 16;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 300 xSize*50 ySize*50])
    set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
    clear xSize sLeft ySize yTop
    %        set(fn,    'windowbuttondownfcn', 'irf_minvar_gui(''ax'')');zoom off;
    ud.h(1)=axes('position',[0.1 0.3 0.5 0.3]); % [x y dx dy]
    ud.h(2)=axes('position',[0.1 0.67 0.5 0.3]); % [x y dx dy]
    ud.h(3)=axes('position',[0.1 0.0 0.5 0.13]); % [x y dx dy]
    linkaxes(ud.h(1:2),'x');
    axis(ud.h(3),'off');
    ud.ht=text(0,1,'','parent',ud.h(3));
    set(fn,'userdata',ud);
    %% initialize probe menu
    hp = uipanel('Title','Probe','FontSize',12,'BackgroundColor',[1 0.95 1],'Position',[.7 .0 .3 .39]);
    inp.probe.type                       = uicontrol('Parent',hp,'String','spherical probe|cylindrical probe|specify probe area','style','popup','Position',[2 230 150 30],'backgroundcolor','white','Callback', @setprobetype);
    surf=lp.photocurrent;probtxt='probe surface';for ii=1:numel(surf),probtxt(end+1:end+1+numel(surf{ii}))=['|' surf{ii}];end
    inp.probe.surface                    = uicontrol('Parent',hp,'String',probtxt,                 'Position',[2 210 150 30],'style','popup','backgroundcolor','white');
    inp.probe.total_vs_sunlit_area_text  = uicontrol('Parent',hp,'String','total/sunlit area',     'Position',[0   185 120 30]);
    inp.probe.total_vs_sunlit_area_value = uicontrol('Parent',hp,'String',num2str(ud.probe.total_vs_sunlit_area),'style','text',...
      'Position',[120 185 70 30],'backgroundcolor','white');
    inp.probe.length_text                = uicontrol('Parent',hp,'String','probe length [cm]',     'Position',[0   160 120 30]);
    inp.probe.length_value               = uicontrol('Parent',hp,'String','','style','text',       'Position',[120 160 70 30],'backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
    inp.probe.radius_text                = uicontrol('Parent',hp,'String','probe radius [cm]',     'Position',[0   135 120 30]);
    inp.probe.radius_value               = uicontrol('Parent',hp,'String',num2str(ud.probe.radius),'Position',[120 135 70 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
    inp.probe.bias_current_text          = uicontrol('Parent',hp,'String','bias current [uA]',     'Position',[0   110 120 30]);
    inp.probe.bias_current_value         = uicontrol('Parent',hp,'String','0','style','edit',      'Position',[120 110 70 30],'backgroundcolor','white','Callback','lp.sweep_gui(''update'')');

    inp.UV_factor_text                   = uicontrol('Parent',hp,'String','UV factor',             'Position',[0   80 60 30]);
    inp.UV_factor_value                  = uicontrol('Parent',hp,'String',num2str(ud.UV_factor),   'Position',[70  80 100 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
    inp.Rsun_text                        = uicontrol('Parent',hp,'String','Rsun [AU]',             'Position',[0   55 60 30]);
    inp.Rsun_value                       = uicontrol('Parent',hp,'String',num2str(ud.R_sun),       'Position',[70  55 100 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
    inp.U_text                           = uicontrol('Parent',hp,'String','U [V]',                 'Position',[0   30 60 30]);
    inp.U_value                          = uicontrol('Parent',hp,'String',ud.U_string,             'Position',[70  30 100 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
    inp.update                           = uicontrol('Parent',hp,'String','Update',                'Position',[0   0 60 30],'Callback','lp.sweep_gui(''update'')');
    inp.reset                            = uicontrol('Parent',hp,'String','Reset',                 'Position',[70  0 60 30],'callback','lp.sweep_gui(''initialize'')');
    %% initialize s/c menu
    hsc = uipanel('Title','Spacecraft','FontSize',12,'BackgroundColor',[.95 1 1],'Position',[.7 .39 .3 .35]);
    inp.flag_sc                                    = uicontrol('Parent',hsc,'style','radio','String','Model spacecraft','Value',0,             'Position',[0   205 120 25]);
    inp.sc.example                                 = uicontrol('Parent',hsc,'String','Example spacecraft|Cluster|Solar Orbiter|THEMIS|Cassini','Position',[0   180 150 25],'style','popup','backgroundcolor','white','Callback', @setscexample);
    surf=lp.photocurrent;probtxt='spacecraft surface';for ii=1:numel(surf),probtxt(end+1:end+1+numel(surf{ii}))=['|' surf{ii}];end
    inp.sc.surface                                 = uicontrol('Parent',hsc,'String',probtxt,                                                  'Position',[0   155 150 30],'style','popup','backgroundcolor','white');
    inp.sc.total_area_text                         = uicontrol('Parent',hsc,'String','Total area [m2]',                                        'Position',[0   135 120 25]);
    inp.sc.total_area_value                        = uicontrol('Parent',hsc,'String',num2str(ud.sc.total_area),'style','edit',                 'Position',[120 135 50 25],'backgroundcolor','white');
    inp.sc.sunlit_area_text                        = uicontrol('Parent',hsc,'String','Sunlit area [m2]',                                       'Position',[0   110 120 25]);
    inp.sc.sunlit_area_value                       = uicontrol('Parent',hsc,'String',num2str(ud.sc.sunlit_area),'style','edit',                'Position',[120 110 50 25],'backgroundcolor','white');
    inp.sc.antenna_guard_area_text                 = uicontrol('Parent',hsc,'String','Sunlit guard area [m2]',                                 'Position',[0   85 120 25],'Tooltipstring','Cross section area of pucks and guards, assuming similar photoelectron emission as antenna');
    inp.sc.antenna_guard_area_value                = uicontrol('Parent',hsc,'String','0','style','edit',                                       'Position',[120 85 50 25],'backgroundcolor','white');
    inp.sc.probe_refpot_as_fraction_of_scpot_text  = uicontrol('Parent',hsc,'String','Probe refpot/scpot',                                     'Position',[0   60 120 25],'Tooltipstring','The ratio between the probe reference potential and satellite potential');
    inp.sc.probe_refpot_as_fraction_of_scpot_value = uicontrol('Parent',hsc,'String',num2str(ud.sc.probe_refpot_as_fraction_of_scpot),         'Position',[120 60 50 25],'style','edit','backgroundcolor','white');
    inp.sc.number_of_probes_text                   = uicontrol('Parent',hsc,'String','Number of probes',                                       'Position',[0   35 120 25]);
    inp.sc.number_of_probes_value                  = uicontrol('Parent',hsc,'String',num2str(ud.sc.number_of_probes),                          'Position',[120 35 50 25],'style','edit','backgroundcolor','white');
    inp.sc.probe_distance_to_spacecraft_text       = uicontrol('Parent',hsc,'String','distance probe-sc [m]',                                  'Position',[0   10 120 25]);
    inp.sc.probe_distance_to_spacecraft_value      = uicontrol('Parent',hsc,'String',num2str(ud.sc.probe_distance_to_spacecraft),              'Position',[120 10 50 25],'style','edit','backgroundcolor','white');
    %% initialize plasma menu
    hpl= uipanel('Title','Plasma','FontSize',12,'BackgroundColor',[1 1 .95],'Position',[.7 .74 .3 .2]);
    inp.n         = uicontrol('Parent',hpl,'String','Ne [cc]','Position',[0 0 80 25]);
    inp.n_value   = uicontrol('Parent',hpl,'String','1','style','edit','Position',[80 0 90 25],'backgroundcolor','white');
    inp.T         = uicontrol('Parent',hpl,'String','T [eV]','Position',[0 25 80 25]);
    inp.T_value   = uicontrol('Parent',hpl,'String','1 1','style','edit','Position',[80 25 90 25],'backgroundcolor','white');
    inp.m         = uicontrol('Parent',hpl,'String','m [mp],0=me','Position',[0 50 80 25]);
    inp.m_value   = uicontrol('Parent',hpl,'String','0 1','style','edit','Position',[80 50 90 25],'backgroundcolor','white');
    inp.q         = uicontrol('Parent',hpl,'String','q [e]','Position',[0 75 80 25]);
    inp.q_value   = uicontrol('Parent',hpl,'String','-1 1','style','edit','Position',[80 75 90 25],'backgroundcolor','white');
    inp.vsc       = uicontrol('Parent',hpl,'String','Vsc [m/s]','Position',[0 100 80 25]);
    inp.vsc_value = uicontrol('Parent',hpl,'String','0','style','edit','Position',[80 100 90 25],'backgroundcolor','white');
    %% initialize plot menu
    hpl= uipanel('Title','Top panel','FontSize',12,'BackgroundColor',[1 1 .95],'Position',[.7 .94 .3 .06]);
    inp.toppanel.plot = uicontrol('Parent',hpl,'String','Resistance|Satellite IU|Antenna noise','Position',[0 0 150 25],'style','popup','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');

    ud.inp=inp;
    set(gcf,'userdata',ud);
    %
    lp.sweep_gui('update');
  case 'update'
    disp('Updating...');
    ud=get(gcf,'userdata');
    %% get input parameters
    inp=ud.inp;
    ud.U=eval(get(inp.U_value,'string'));
    ud.UV_factor=str2double(get(inp.UV_factor_value,'string'));
    ud.R_sun=str2double(get(inp.Rsun_value,'string'));
    new_radius=str2double(get(inp.probe.radius_value,'string'));
    if (new_radius ~= ud.probe.radius)
      ud.probe.radius=new_radius;
      ud.probe.capacitance=0; % force recalculate capacitanc
    end
    new_length=str2double(get(inp.probe.length_value,'string'));
    if (new_length ~= ud.probe.length)
      ud.probe.length=new_length;
      ud.probe.capacitance=0; % force recalculate capacitanc
    end
    ud.probe.bias_current=str2double(get(inp.probe.bias_current_value,'string'))*1e-6; % convert from uA to A
    surf=lp.photocurrent;
    ud.probe.surface=surf{max(1,get(inp.probe.surface,'Value')-1)}; %
    ud.n  =eval(['[' get(inp.n_value,'string')   ']' ]);
    ud.T  =eval(['[' get(inp.T_value,'string')   ']' ]);
    ud.m  =eval(['[' get(inp.m_value,'string')   ']' ]);
    ud.q  =eval(['[' get(inp.q_value,'string')   ']' ]);
    ud.vsc=eval(['[' get(inp.vsc_value,'string') ']' ]);
    % if scflag read in sc parameters
    ud.flag_use_sc=get(inp.flag_sc,'Value');
    ud.toppanel=get(inp.toppanel.plot,'Value');
    if ud.flag_use_sc
      ud.probe_refpot_as_fraction_of_scpot=str2double(get(inp.sc.probe_refpot_as_fraction_of_scpot_value,'string'));
      ud.sc.number_of_probes=str2double(get(inp.sc.number_of_probes_value,'string')); % in cm
      ud.sc.cross_section_area=str2double(get(inp.sc.sunlit_area_value,'string'));
      ud.sc.total_area=str2double(get(inp.sc.total_area_value,'string'));
      ud.sc.capacitance=4*pi*Units.eps0*sqrt(ud.sc.total_area/4/pi); % assume sphere having the specified total area
      ud.sc.type='spherical';
      surf=lp.photocurrent;
      ud.sc.surface=surf{max(1,get(inp.sc.surface,'Value')-1)}; %
      ud.sc.antenna_guard_area=str2double(get(inp.sc.antenna_guard_area_value,'string'));
    end
    %% calculate IU curves
    Upot=ud.U(:);
    probe=ud.probe;
    switch ud.probe.type
      case 'spherical'
        probe.cross_section_area=pi*(ud.probe.radius*.01)^2;
        probe.total_area=4*probe.cross_section_area;
        if ~isfield(probe,'capacitance') || probe.capacitance==0 % if not defined or put to zero, calculate
          probe.capacitance=4*pi*Units.eps0*ud.probe.radius*.01;
        end
      case 'cylindrical'
        probe.cross_section_area=2*ud.probe.radius*ud.probe.length*0.01^2;
        probe.total_area=pi*probe.cross_section_area;
        probe.length=str2double(get(inp.probe.length_value,'string'));
        if ~isfield(probe,'capacitance') || probe.capacitance==0 % if not defined or put to zero, calculate
          %                    probe.capacitance=2*pi*Units.eps0*ud.probe.length*0.01/log(ud.probe.length/ud.probe.radius); % assuming length >> radius
          probe.capacitance=irf_estimate('capacitance_cylinder',ud.probe.radius*0.01,ud.probe.length*0.01); % assuming length >> radius
        end
      case 'arbitrary'
        probe.cross_section_area=str2double(get(inp.probe.radius_value,'string'))*0.01^2;
        probe.total_area=str2double(get(inp.probe.length_value,'string'))*0.01^2;
        set(inp.probe.total_vs_sunlit_area_value,'string',num2str(probe.total_area/probe.cross_section_area,3));
        if ~isfield(probe,'capacitance') || probe.capacitance==0 % if not defined or put to zero, calculate
          probe.capacitance=4*pi*Units.eps0*sqrt(probe.total_area/4/pi); % assume sphere having the specified total area
        end
    end
    ud.probe=probe;
    [J_probe,J_photo,J_plasma]=lp.probe_current(probe,Upot,ud.R_sun,ud.UV_factor,ud);
    dUdI=gradient(Upot,J_probe);
    ud.I=J_probe;
    ud.dUdI=dUdI;
    % if scflag then calculate s/c IU curve
    if ud.flag_use_sc
      J_sc=lp.probe_current(ud.sc,Upot,ud.R_sun,ud.UV_factor,ud);
      ud.I_sc=J_sc;
      ud.U_sc=Upot;
      ud.dUdI_sc=gradient(Upot,J_sc);
    end
    % if scflag calculate probe to sc IU curve
    if ud.flag_use_sc
      % reduce probe curve to reasonable number of points (derivative does
      % not change more than 10% between points
      ind=ones(size(Upot));
      ilast=1;
      slopechange=0.2; % how long slope change is allowed
      dudi=gradient(Upot,J_probe);
      for ii=2:numel(Upot)
        if abs((dudi(ii)-dudi(ilast))/dudi(ilast))<slopechange && Upot(ii)-Upot(ilast)<.4
          ind(ii)=NaN;
        else % keep previous point as last
          ind(ii-1)=1;
          ilast=ii-1;
        end
      end
      Ibias=J_probe(ind==1); % bias current to probe measured with respect to sc
      Ubias=Upot(ind==1);

      %
      antena_guard_area_factor=(ud.probe.cross_section_area+ud.sc.antenna_guard_area/ud.sc.number_of_probes)/ud.probe.cross_section_area;
      Iprobe=Ibias*antena_guard_area_factor;
      % zero approximation
      Isat=ud.sc.number_of_probes*Iprobe;
      Usatsweep=interp1(J_sc,Upot,Isat/5,'nearest'); % floating potential of sc during sweep, assuming of only 1/5 of bias current electrons esacpe to space
      Isat_probe_photoelectrons=Isat*0; % electrons from probes hitting spacecraft
      Uprobe2plasma=zeros(size(Ibias)); % initialize
      Uproberefsweep=zeros(size(Ibias)); % initialize
      Jprobephotoreturn=zeros(size(Ibias)); % initialize
      FitError=zeros(size(Ibias)); % initialize
      refpotvec=Ubias(1:4:end)*ud.probe_refpot_as_fraction_of_scpot;
      [probepotgrid,refpotgrid] = meshgrid(Ubias,refpotvec);
      Jprobephotogrid=lp.probe_current(probe,probepotgrid-refpotgrid,ud.R_sun,ud.UV_factor,[]);
      Jprobephotoescapingscgrid=lp.probe_current(probe,probepotgrid,ud.R_sun,ud.UV_factor,[]);
      Jprobeplasmagrid=lp.probe_current(probe,probepotgrid,ud.R_sun,0,ud);
      Jprobegrid=Jprobeplasmagrid+Jprobephotogrid;
      Jprobephotoe2scgrid=Jprobephotogrid-Jprobephotoescapingscgrid;
      for ii=1:numel(Iprobe)
        % plasma current with UV factor zero
        satpot=Usatsweep(ii);
        uprobe=interp1(J_probe,Upot,Ibias(ii))+satpot*ud.probe_refpot_as_fraction_of_scpot;
        ibias=Ibias(ii);
        jj=1;
        err=1;
        while err > .001 && jj<10
          refpot=satpot*ud.probe_refpot_as_fraction_of_scpot;
          if isnan(refpot), break; end
          uprobe=interp1(interp1(refpotvec,Jprobegrid,refpot,'linear','extrap'),Ubias,ibias);
          if isnan(uprobe); uprobe=interp1(J_probe,Upot,Ibias(ii))+satpot*ud.probe_refpot_as_fraction_of_scpot;end
          J_probe_photoe2sc=interp2(probepotgrid,refpotgrid,Jprobephotoe2scgrid,uprobe,refpot);
          new_Isat_probe_photoelectrons=-J_probe_photoe2sc*antena_guard_area_factor*ud.sc.number_of_probes;
          err=abs(new_Isat_probe_photoelectrons-Isat_probe_photoelectrons(ii))/abs(new_Isat_probe_photoelectrons);
          if isnan(err), err=0; end
          Isat_probe_photoelectrons(ii)=new_Isat_probe_photoelectrons;
          satpot=interp1(J_sc,Upot,-Isat(ii)-Isat_probe_photoelectrons(ii));
          jj=jj+1;
        end
        FitError(ii)=err;
        Jprobephotoreturn(ii)=J_probe_photoe2sc;
        Usatsweep(ii)=satpot;
        Uprobe2plasma(ii)=uprobe;
        Uproberefsweep(ii)=refpot;
      end

      Uprobe2sc        =Uprobe2plasma-Usatsweep;
      Uprobe2refpot    =Uprobe2plasma-Uproberefsweep;
      dUdI_probe2plasma=gradient(Uprobe2plasma,Ibias);
      dUdI_probe2sc    =gradient(Uprobe2sc,Ibias);
      dUdI_probe2refpot=gradient(Uprobe2refpot,Ibias);
      dUdI             =dUdI_probe2refpot;
      Upot             =Uprobe2refpot;
      J_probe          =Ibias;
    end
    %% plot IU curve
    info_txt='';
    h=ud.h;
    flag_add_bias_point_values=0; % default
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BOTTOM PANEL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(h(1),Upot,J_probe*1e6,'k');
    grid(h(1),'on');
    xlabel(h(1),'U [V]');
    ylabel(h(1),'I [\mu A]');
    hold(h(1),'on');
    if ud.flag_use_sc % add probe potential wrt s/c and plasma
      plot(h(1),Ubias,Ibias*1e6,'k','linewidth',0.5);
      plot(h(1),Uprobe2plasma,Ibias*1e6,'r.','linewidth',1.5);
      plot(h(1),Uprobe2sc,Ibias*1e6,'b','linewidth',1.5);
      plot(h(1),Uprobe2sc,Jprobephotoreturn*1e6,'color',[0 0.5 0],'linewidth',1.5);
      plot(h(1),Usatsweep,Ibias*1e6,'color',[0.5 0 0.5],'linewidth',1);
      irf_legend(h(1),'probe to plasma',[1 1.01],'color','r');
      irf_legend(h(1),'probe to reference',[1 1.05],'color','k');
      irf_legend(h(1),'probe to spacecraft (bias)',[0.02 1.01],'color','b');
      irf_legend(h(1),'probe photo e- to s/c',[0.02 1.05],'color',[0 0.5 0]);
      irf_legend(h(1),'Satellite potential',[0.98 0.03],'color',[0.5 0 0.5]);
      axis(h(1),'auto x');
      if ud.probe.bias_current ~= 0 && -ud.probe.bias_current>min(Ibias) && -ud.probe.bias_current<max(Ibias) % draw bias current line
        plot(h(1),[Uprobe2sc(1) Uprobe2plasma(end)],ud.probe.bias_current.*[-1 -1].*1e6,'k-.','linewidth',0.5);
        text(Uprobe2sc(1),ud.probe.bias_current*(-1)*1e6,'bias','parent',h(1),'horizontalalignment','left','verticalalignment','bottom');
        flag_add_bias_point_values=1;
      end
    else % add photoelectron and photoelectron currents
      plot(h(1),Upot,J_photo*1e6,'r','linewidth',0.5);
      irf_legend(h(1),'    total',      [0.98 0.03],'color','k');
      irf_legend(h(1),' photoelectrons',[0.98 0.08],'color','r');
      clr=[0.5 0 0; 0 0.5 0; 0 0 0.5];
      for ii=1:length(J_plasma)
        plot(h(1),Upot,J_plasma{ii}*1e6,'linewidth',.5,'color',clr(:,ii));
        irf_legend(h(1),['plasma ' num2str(ii)],[0.98 0.08+ii*0.05],'color',clr(:,ii));
      end
      if ud.probe.bias_current ~= 0 && -ud.probe.bias_current>min(J_probe) && -ud.probe.bias_current<max(J_probe) % draw bias current line
        plot(h(1),[Upot(1) Upot(end)],ud.probe.bias_current.*[-1 -1].*1e6,'k-.','linewidth',0.5);
        text(Upot(1),ud.probe.bias_current*(-1)*1e6,'-bias','parent',h(1),'horizontalalignment','left','verticalalignment','bottom');
        flag_add_bias_point_values=1;
      end
    end
    hold(h(1),'off');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFORMATION PANEL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Rmin = min(abs(dUdI)); % minimum resistance
    fcr=1/2/pi/Rmin/probe.capacitance;
    disp(['Rmin=' num2str(Rmin,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, f_{CR}=' num2str(fcr,3) 'Hz.']);
    if ud.flag_use_sc
      info_txt=[info_txt '\newline probe to plasma Rmin=' num2str(min(abs(dUdI_probe2plasma)),3) ' Ohm'];
      info_txt=[info_txt '\newline probe to reference Rmin=' num2str(min(abs(dUdI)),3) ' Ohm'];
      info_txt=[info_txt '\newline probe to spacecraft Rmin=' num2str(min(abs(dUdI_probe2sc)),3) ' Ohm'];
    else
      info_txt=[info_txt '\newline Rmin=' num2str(Rmin,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.'];
    end

    if min(J_probe)<0 && max(J_probe)>0                   % display information on Ufloat
      Ufloat=interp1(J_probe,Upot,0); % floating potential
      ii=isfinite(Upot);Rfloat=interp1(Upot(ii),dUdI(ii),Ufloat);
      fcr=1/2/pi/Rfloat/probe.capacitance;
      info_txt=[info_txt '\newline Probe: Ufloat=' num2str(Ufloat,3) 'V, Rfloat= ' num2str(Rfloat,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.'];
      disp(['Probe: Ufloat=' num2str(Ufloat,3) ' V, Rfloat=' num2str(Rfloat,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.']);
    end
    if flag_add_bias_point_values
      Ubias=interp1(J_probe,Upot,-ud.probe.bias_current); % floating potential
      ii=isfinite(Upot);
      Rbias=interp1(Upot(ii),dUdI(ii),Ubias);
      fcr=1/2/pi/Rbias/probe.capacitance;
      disp(['Rbias=' num2str(Rbias,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.']);
      info_txt=[info_txt '\newline Probe: (without s/c) Ubias=' num2str(Ubias,3)  ', Rbias=' num2str(Rbias,3) 'Ohm, fcr=' num2str(fcr,3) 'Hz.'];
      if ud.flag_use_sc
        Uscbias=interp1(J_probe,Usatsweep,-ud.probe.bias_current); % floating potential
        ii=isfinite(Upot);
        Rscbias=interp1(ud.U_sc(ii),ud.dUdI_sc(ii),Uscbias);
        fcr=1/2/pi/Rscbias/ud.sc.capacitance;
        disp(['Spacecraft (biased): Rbias=' num2str(Rbias,3) ' Ohm, C=' num2str(ud.sc.capacitance*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.']);
        info_txt=[info_txt '\newline Spacecraft (biased): Ubias=' num2str(Uscbias,3)  ', Rbias=' num2str(Rscbias,3) 'Ohm, fcr=' num2str(fcr,3) 'Hz.'];
        Usp=-interp1(Ibias,Uprobe2sc,-ud.probe.bias_current); % floating potential
        disp(['Spacecraft to Probe: Usp=' num2str(Usp,3) 'V.']);
        info_txt=[info_txt '\newline Spacecraft to Probe: Usp=' num2str(Usp,3)  ' V.'];
      end
    end
    if ud.flag_use_sc && min(ud.I_sc)<0 && max(ud.I_sc)>0 % display information on Ufloat
      Uscfloat=interp1(ud.I_sc,ud.U_sc,0); % floating potential
      ii=isfinite(ud.U_sc);Rscfloat=interp1(ud.U_sc(ii),ud.dUdI_sc(ii),Uscfloat);
      info_txt=[info_txt '\newline Spacecraft: Ufloat=' num2str(Uscfloat,3) 'V, Rfloat= ' num2str(Rscfloat,3) ' Ohm, C=' num2str(ud.sc.capacitance*1e12,3) 'pF'];
      disp(['Spacecraft: Ufloat=' num2str(Uscfloat,3) ' V, Rfloat=' num2str(Rscfloat,3) ' Ohm, C=' num2str(ud.sc.capacitance*1e12,3) 'pF']);
    end
    if ud.UV_factor>0                                     % display photoelectron saturation current
      info_txt=[info_txt '\newline Probe: photo e- Io = ' num2str(ud.UV_factor*lp.photocurrent(1,-1,ud.R_sun,ud.probe.surface)*1e6,3) '[\mu A/m^2]'];
      if ud.R_sun~=1
        info_txt=[info_txt '  (' num2str(ud.UV_factor*lp.photocurrent(1,-1,1,ud.probe.surface)*1e6,3) ' \mu A/m^2 at 1 AU)'];
      end
      if ud.flag_use_sc
        info_txt=[info_txt '\newline Spacecraft: photo e- Io = ' num2str(ud.UV_factor*lp.photocurrent(1,-1,ud.R_sun,ud.sc.surface)*1e6,3) '[\mu A/m^2]'];
        if ud.R_sun~=1
          info_txt=[info_txt '  (' num2str(ud.UV_factor*lp.photocurrent(1,-1,1,ud.sc.surface)*1e6,3) ' \mu A/m^2 at 1 AU)'];
        end
      end
    end

    axis(h(3),'off');
    set(ud.ht,'string',info_txt,'fontsize',11);
    set(gcf,'userdata',ud);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TOP PANEL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ud.toppanel==1 % plot resistance
      plot(h(2),Upot,dUdI,'k');
      grid(h(2),'on');xlabel(h(2),'U [V]');
      ylabel(h(2),'dU/dI [\Omega]');
      if ud.flag_use_sc % add probe resistance wrt plasma and s/c
        hold(h(2),'on');
        plot(h(2),Uprobe2plasma,dUdI_probe2plasma,'r','linewidth',1.5);
        plot(h(2),Uprobe2sc,dUdI_probe2sc,'b','linewidth',1.5);
        hold(h(2),'off');
      end
      axis(h(2),'auto y');
      set(h(2),'yscale','log')
      linkaxes(ud.h(1:2),'x');
    elseif ud.toppanel==2 && ud.flag_use_sc % plot spacecraft IU
      plot(h(2),ud.U_sc,ud.I_sc*1e6,'k');
      grid(h(2),'on');xlabel(h(2),'U [V]');
      ylabel(h(2),'I [\mu A/m^2]');
      linkaxes(ud.h(1:2),'x');
    elseif ud.toppanel==3 % plot probe noise levels
      E_int_range=[1e-18 1e-8]; % range of E intensity
      f_range=[1e-3 9.9e5];      % frequencies in Hz
      f_noise_range=10.^(log10(f_range(1)):.1:log10(f_range(2)));
      f=f_noise_range';
      ud.sc.probe_distance_to_spacecraft=str2double(get(inp.sc.probe_distance_to_spacecraft_value,'string'));
      antenna_eff_length=ud.sc.probe_distance_to_spacecraft*2; % efficient length of antenna from satellite center
      C_antenna=ud.probe.capacitance;
      n=ud.n*1e6;                                             % TEMPORAR SOLUTION!!!!
      T_plasma=ud.T(1);                                       % TEMPORAR SOLUTION!!!!
      A_antenna=ud.probe.total_area;
      distance_to_Sun=ud.R_sun; % AU
      UV=ud.UV_factor;

      noise_preamp=10e-9; % preamplifier noise 10nV/Hz1/2
      noise_preamp_level=(noise_preamp/antenna_eff_length)^2;
      f_break=400; % transition frequency at which 1/f noise is starting

      % instrumental noise
      noise_instr_X=[f_range(1) f_break f_range(2)];
      noise_instr_Y=noise_preamp_level*[sqrt(f_break/f_range(1)) 1  1];

      if exist('Rbias','var')
        T_eV=1; % photoelectron temperature
        noise_thermal_bias=4*Units.e*T_eV*sqrt(Rbias^2./(1+(2*pi*f).^2*Rbias^2*C_antenna^2))/antenna_eff_length^2;
        nu=n/2*sqrt(8*Units.kB*T_plasma/pi/Units.me)*A_antenna;
        noise_shot_plasma_bias=2*nu*Units.e^2*Rbias^2./(1+(2*pi*f).^2*Rbias^2*C_antenna^2)/antenna_eff_length^2;
        Iphoto=abs(lp.probe_current(probe,-1,distance_to_Sun,UV,[]));
        nu=Iphoto/Units.e;
        noise_shot_photoelectrons_bias=2*nu*Units.e^2*Rbias^2./(1+(2*pi*f).^2*Rbias^2*C_antenna^2)/antenna_eff_length^2;
        noise_total_bias=noise_shot_photoelectrons_bias+noise_shot_plasma_bias+noise_thermal_bias;
      end
      if exist('Rfloat','var')
        T_eV=1; % photoelectron temperature
        noise_thermal_nobias=4*Units.e*T_eV*sqrt(Rfloat^2./(1+(2*pi*f).^2*Rfloat^2*C_antenna^2))/antenna_eff_length^2;
        nu=n/2*sqrt(8*Units.kB*T_plasma/pi/Units.me)*A_antenna;
        noise_shot_plasma_nobias=2*nu*Units.e^2*Rfloat^2./(1+(2*pi*f).^2*Rfloat^2*C_antenna^2)/antenna_eff_length^2;
        Iphoto=abs(lp.probe_current(probe,-1,distance_to_Sun,UV,[]));
        nu=Iphoto/Units.e;
        noise_shot_photoelectrons_nobias=2*nu*Units.e^2*Rfloat^2./(1+(2*pi*f).^2*Rfloat^2*C_antenna^2)/antenna_eff_length^2;
        noise_total_nobias=noise_shot_photoelectrons_nobias+noise_shot_plasma_nobias+noise_thermal_nobias;
      end
      if 1 % plot
        linkaxes(ud.h(1:2),'off');hold(h(2),'off');
        plot(h(2),noise_instr_X,noise_instr_Y,'k');
        set(h(2),'xscale','log','yscale','log','xlim',f_range,'ylim',E_int_range);
        set(h(2),'xtick',10.^( floor(log10(f_range(1))):ceil(log10(f_range(2))) ) );
        set(h(2),'ytick',10.^( floor(log10(E_int_range(1))):ceil(log10(E_int_range(2))) ) );
        set(h(2),'MinorGridLineStyle','none','FontSize',12)
        grid(h(2),'on');
        ht=text(1e3,noise_preamp_level*1.5,'instrument');set(ht,'fontsize',14,'color','k');
        hold(h(2),'on');
        ax=h(2);
        if exist('Rfloat','var')
          plot(h(2),f,noise_shot_plasma_nobias,'k:')
          plot(h(2),f,noise_shot_photoelectrons_nobias,'LineStyle',':','Color',[0 .5 0])
          plot(h(2),f,noise_thermal_nobias,'b:', ...
            f,noise_total_nobias,'r:')
          ht=text(2e2,noise_thermal_nobias(1)*1.5,'plasma thermal noise (nobias)','parent',ax);set(ht,'fontsize',14,'color','b');
          ht=text(f_range(1),noise_shot_plasma_nobias(1),'shot noise plasma (nobias)','parent',ax);
          set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','left','color','k');
          ht=text(f_range(1),noise_shot_photoelectrons_nobias(1),'shot noise photoelectrons (nobias)','parent',ax);
          set(ht,'fontsize',14,'verticalalignment','top','horizontalalignment','left','color',[0 0.5 0]);
          ht=text(f_range(1),noise_total_nobias(1),'total noise (nobias)','parent',ax);
          set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','left','color','r');
        end
        if exist('Rbias','var')
          plot(h(2),f,noise_shot_plasma_bias,'k', ...
            f,noise_shot_photoelectrons_bias,'g', ...
            f,noise_thermal_bias,'b',...
            f,noise_total_bias,'r');
          ht=text(2e2,noise_thermal_bias(1)*1.5,'plasma thermal noise (bias)','parent',ax);set(ht,'fontsize',14,'color','b');
          ht=text(f_range(1),noise_shot_plasma_bias(1),'shot noise plasma (bias)','parent',ax);
          set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','left');
          ht=text(f_range(1),noise_shot_photoelectrons_bias(1),'shot noise photoelectrons (bias)','parent',ax);
          set(ht,'fontsize',14,'verticalalignment','top','horizontalalignment','left','color',[0 0.5 0]);
          ht=text(f_range(1),noise_total_bias(1),'total noise (bias)','parent',ax);
          set(ht,'fontsize',14,'verticalalignment','bottom','horizontalalignment','left','color','r');
        end

        ylabel(h(2),'Electric field intensity [V^2/m^2/Hz]');
        %ylabel('Electric field intensity [dB V/m/Hz^{1/2}]');
        xlabel(h(2),'Frequency [Hz]');
        hold(h(2),'off');
      end
    end

end
disp('READY!')
end

function setprobetype(hObj,event) %#ok<INUSD>
val = get(hObj,'Value');
data=get(gcf,'userdata');
set(data.inp.probe.radius_text,'string','probe radius [cm]');
set(data.inp.probe.length_text,'string','probe length [cm]');
if val ==1
  data.probe.type='spherical';
  set(data.inp.probe.length_value,'style','text','string','');
elseif val == 2
  data.probe.type='cylindrical';
  set(data.inp.probe.length_value,'style','edit','string',num2str(data.probe.length));
elseif val == 3
  set(data.inp.probe.radius_text,'string','Sunlit area [cm2]');
  set(data.inp.probe.length_text,'string','Total area [cm2]');
  set(data.inp.probe.length_value,'string','1','style','edit');
  set(data.inp.probe.radius_value,'string','.2','style','edit');
  data.probe.type='arbitrary';
end
set(gcf,'userdata',data);
end
function setscexample(hObj,event) %#ok<INUSD>
Units=irf_units;
val = get(hObj,'Value');
data=get(gcf,'userdata');
if val ==1 % do nothing, shows in menu 'Example spacecraft'
elseif val ==2 % Cluster s/c diameter 2.9, height 1.3m
  data.probe.type='spherical';
  set(data.inp.probe.type,'Value',1);
  data.probe.surface='themis';
  set(data.inp.probe.surface,'Value',find(strcmp('cluster',lp.photocurrent))+1);
  set(data.inp.probe.length_value,'style','text','string','');
  set(data.inp.probe.radius_value,'string','4');
  set(data.inp.sc.surface,'Value',find(strcmp('solar cells',lp.photocurrent))+1); % solar cells
  set(data.inp.sc.total_area_value,'string','25.66');
  set(data.inp.sc.sunlit_area_value,'string','3.87');
  set(data.inp.sc.antenna_guard_area_value,'string','0.039');
  set(data.inp.sc.probe_refpot_as_fraction_of_scpot_value,'string','.2');
  set(data.inp.sc.number_of_probes_value,'string','4');
  set(data.inp.sc.probe_distance_to_spacecraft_value,'string','44');
  set(data.inp.Rsun_value,'string','1');
  set(data.inp.probe.total_vs_sunlit_area_value,'string','4');
  data.probe.total_vs_sunlit_area=4;
  set(data.inp.n_value,'string','1');
  set(data.inp.T_value,'string','100 500');
elseif val == 3 % Solar Orbiter with solar panel back side
  data.probe=solar_orbiter('probe');
  data.sc=solar_orbiter('sc');
  p=solar_orbiter('plasma');
  data.plasma=p.perihelion;
  set(data.inp.probe.type,'Value',2); % cylindrical
  set(data.inp.probe.surface,'Value',find(strcmp(data.probe.surface,lp.photocurrent))+1);
  set(data.inp.probe.radius_value,'string',num2str(data.probe.radius/Units.cm,4)); %cm
  set(data.inp.probe.length_value,'style','edit','string',num2str(data.probe.length/Units.cm,4)); % cm
  set(data.inp.probe.total_vs_sunlit_area_value,'string',num2str(data.probe.total_vs_sunlit_area,4));
  set(data.inp.sc.sunlit_area_value,'string',num2str(data.sc.cross_section_area,4));
  set(data.inp.sc.antenna_guard_area_value,'string',num2str(data.sc.antenna_guard_area,4));
  set(data.inp.sc.total_area_value,'string',num2str(data.sc.total_area,4));
  set(data.inp.sc.probe_refpot_as_fraction_of_scpot_value,'string',num2str(data.sc.probe_refpot_as_fraction_of_scpot,4));
  set(data.inp.sc.number_of_probes_value,'string',num2str(data.sc.number_of_probes,4));
  set(data.inp.Rsun_value,'string','0.28');
  set(data.inp.q_value,'string',num2str(data.plasma.q,'%5.0f '));
  set(data.inp.vsc_value,'string',num2str(data.plasma.vsc,'%5.0e '));
  set(data.inp.n_value,'string',num2str(data.plasma.n,'%5.0f '));
  set(data.inp.T_value,'string',num2str(data.plasma.T,'%5.0f '));
elseif val == 4 % THEMIS
  data.probe.type='spherical';
  set(data.inp.probe.type,'Value',1);
  data.probe.surface='themis';
  set(data.inp.probe.surface,'Value',2);
  set(data.inp.probe.length_value,'style','text','string','');
  set(data.inp.probe.radius_value,'string','4');
  set(data.inp.sc.sunlit_area_value,'string','0.71');
  set(data.inp.sc.antenna_guard_area_value,'string','0.0085');
  set(data.inp.sc.total_area_value,'string','4.2');
  set(data.inp.sc.probe_refpot_as_fraction_of_scpot_value,'string','.2');
  set(data.inp.sc.number_of_probes_value,'string','4');
  set(data.inp.sc.probe_distance_to_spacecraft_value,'string','25');
  set(data.inp.Rsun_value,'string','1');
  set(data.inp.probe.total_vs_sunlit_area_value,'string','4');
  data.probe.total_vs_sunlit_area=4;
  set(data.inp.n_value,'string','1');
  set(data.inp.T_value,'string','100 500');
elseif val == 5 % Cassini, sensor - 50mm sphere on 10.9 cm stub (diameter 6.35 mm) (efficient radius 56.5 taking into account stub)
  % sc - >6.7 metres high and >4 metres wide
  data.probe.type='spherical';
  set(data.inp.probe.type,'Value',1);
  data.probe.surface='cassini';
  set(data.inp.probe.surface,'Value',3);
  set(data.inp.probe.length_value,'style','text','string','');
  set(data.inp.probe.radius_value,'string','2.5');
  set(data.inp.sc.sunlit_area_value,'string','27');
  set(data.inp.sc.antenna_guard_area_value,'string','0.0');
  set(data.inp.sc.total_area_value,'string','140');
  set(data.inp.sc.probe_refpot_as_fraction_of_scpot_value,'string','.2');
  set(data.inp.sc.number_of_probes_value,'string','1');
  set(data.inp.Rsun_value,'string','9.5');
  set(data.inp.probe.total_vs_sunlit_area_value,'string','4');
  data.probe.total_vs_sunlit_area=4;
  set(data.inp.n_value,'string','0.1');
  set(data.inp.T_value,'string','1 1');
end
set(gcf,'userdata',data);
end
function setplasmaexample(hObj,event) %#ok<DEFNU,INUSD>
Units=irf_units; %#ok<NASGU>
val = get(hObj,'Value');
data=get(gcf,'userdata');
if val ==1 % do nothing, shows in menu 'Example spacecraft'
elseif val ==2 % solar wind at 1AU
  set(data.inp.Rsun_value,'string','1');
  set(data.inp.probe.total_vs_sunlit_area_value,'string','4');
  data.probe.total_vs_sunlit_area=4;
  set(data.inp.n_value,'string','1');
  set(data.inp.T_value,'string','100 500');
end
set(gcf,'userdata',data);
end
