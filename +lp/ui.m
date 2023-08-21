classdef ui < handle
  %LP.UI User interface to +lp package
  %   LP.UI(spacecraft,probe,plasma,parameters)
  properties (SetAccess = protected)
    SpacecraftList % Nr 1 is always user defined
    spacecraftUsed
    ProbeList % Nr 1 is always user defined
    probeUsed % 1..N from ProbeList
    ProbeSurfaceList % Nr 1 is always user defined
    ProbeSurfacePhotoemission % in A/m^2 units
    probeSurfaceUsed % 1..N from ProbeSurfaceList
    PlasmaList
    plasmaUsed
    figHandle
    Axes     % top,bottom,infotext
    UserData % user data
    InputParameters =struct('factorUV',1,'rSunAU',1,...
      'vectorUString','-5:0.2:40','vectorU',-5:0.2:40,'biasCurrent',0)
    Output
  end
  methods
    function obj = ui(varargin)
      %% Check input
      args=varargin;
      while ~isempty(args)
        if isa(args{1},'lp.spacecraft')
          obj.SpacecraftList = args{1};
          args(1) =[];
        elseif isa(args{1},'lp.lprobe')
          obj.ProbeList = args{1};
          args(1) =[];
        elseif isa(args{1},'lp.plasma')
          obj.PlasmaList = args{1};
          args(1) =[];
        else
          error('lp.ui input unknown!');
        end
      end
      if ~isa(obj.SpacecraftList,'lp.spacecraft')
        obj.SpacecraftList = lp.default_spacecraft;
      end
      if ~isa(obj.ProbeList,'lp.probe')
        obj.ProbeList = lp.default_lprobe;
      end
      if ~isa(obj.PlasmaList,'lp.plasma')
        obj.PlasmaList = lp.default_plasma;
      end
      obj.SpacecraftList = [obj.SpacecraftList(1) obj.SpacecraftList];
      obj.SpacecraftList(1).name = 'user defined';
      obj.ProbeList = [obj.ProbeList(1) obj.ProbeList];
      obj.ProbeList(1).name = 'user defined';
      obj.PlasmaList = [obj.PlasmaList(1) obj.PlasmaList];
      obj.PlasmaList(1).name = 'user defined';
      obj.ProbeSurfaceList = [{'user defined'} lp.photocurrent];
      obj.ProbeSurfacePhotoemission = zeros(size(obj.ProbeSurfaceList));
      for j = 2:numel(obj.ProbeSurfacePhotoemission)
        obj.ProbeSurfacePhotoemission(j) = lp.photocurrent(1,0,1,obj.ProbeSurfaceList{j});
      end

      %% Initialize IDE
      obj.new_ide();
      obj.spacecraftUsed = 2;
      obj.probeUsed = 2;
      obj.plasmaUsed = 2;
      obj.set_plasma_model(obj.plasmaUsed);
      obj.set_probe_model(  obj.probeUsed);
      obj.get_u_interval;
    end
    function new_ide(obj)
      %% initialize figure
      set(0,'defaultLineLineWidth', 1.5);
      figH=figure;
      obj.figHandle=figH;
      clf reset;
      clear h;
      set(figH,'color','white'); % white background for figures (default is grey)
      set(figH,'PaperUnits','centimeters')
      set(figH,'defaultAxesFontSize',14);
      set(figH,'defaultTextFontSize',14);
      set(figH,'defaultAxesFontUnits','pixels');
      set(figH,'defaultTextFontUnits','pixels');
      xSize = 13; ySize = 16;
      xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
      set(figH,'PaperPosition',[xLeft yTop xSize ySize])
      set(figH,'Position',[100 300 xSize*50 ySize*50])
      set(figH,'paperpositionmode','auto') % to get the same printing as on screen
      clear xSize sLeft ySize yTop
      %        set(fn,    'windowbuttondownfcn', 'irf_minvar_gui(''ax'')');zoom off;
      obj.Axes.bottom = axes('position',[0.1 0.3 0.5 0.3]); % [x y dx dy]
      obj.Axes.top    = axes('position',[0.1 0.67 0.5 0.3]); % [x y dx dy]
      linkaxes([obj.Axes.top obj.Axes.bottom],'x');

      obj.Axes.infotext = axes('position',[0.1 0.0 0.5 0.13]); % [x y dx dy]
      axis(obj.Axes.infotext,'off');
      obj.Axes.infoTextHandle = text(0,1,'','parent',obj.Axes.infotext);
      %% initialize probe menu
      colPanelBg = [1 0.95 1];
      hp         = uipanel('Title','Probe','FontSize',12,'BackgroundColor',colPanelBg,'Position',[.7 .0 .3 .39]);
      uiPar      = {'Parent',hp,'backgroundcolor',colPanelBg,'style','text'};
      popupProbeTxt        = obj.popup_list(obj.ProbeList);
      popupProbeSurfaceTxt = obj.popup_list(obj.ProbeSurfaceList);
      inp.probe.type.text                 = uicontrol(uiPar{:},'String','type','style','text',   'Position',[0   230 60   20]);
      inp.probe.type.value                = uicontrol(uiPar{:},'String',popupProbeTxt,           'Position',[60  230 130  20],'style','popup','Callback',@(src,evt)obj.set_probe_model(src,evt));
      inp.probe.surface.text              = uicontrol(uiPar{:},'String','surface',               'Position',[0   210 60   20]);
      inp.probe.surface.value             = uicontrol(uiPar{:},'String',popupProbeSurfaceTxt,    'Position',[60  210 130  20],'style','popup','Callback',@(src,evt)obj.set_probe_surface(src,evt));
      inp.probe.surfacePhotoemission.text = uicontrol(uiPar{:},'String','Iph @1AU,UV=1 [uA/m2]', 'Position',[0   190 120 20],'style','text');
      inp.probe.surfacePhotoemission.value= uicontrol(uiPar{:},'String','0',                     'Position',[120 190 70  20],'style','edit','Callback',@(src,evt)obj.get_probe_surface);
      inp.probe.areaTotalVsSunlit.text    = uicontrol(uiPar{:},'String','total/sunlit area',     'Position',[0   170 120 20],'style','text');
      inp.probe.areaTotalVsSunlit.value   = uicontrol(uiPar{:},'String','',                      'Position',[120 170 70  20]);
      inp.probe.radiusSphere.text         = uicontrol(uiPar{:},'String','sphere radius [cm]',    'Position',[0   150 120 20]);
      inp.probe.radiusSphere.value        = uicontrol(uiPar{:},'String','','style','edit',       'Position',[120 150 70  20],'backgroundcolor','white','Callback',@(src,evt)obj.get_probe_radius_sphere);
      inp.probe.radiusSphere.SIconversion = 1e-2;
      inp.probe.lengthWire.text           = uicontrol(uiPar{:},'String','cyl/wire length [cm]',  'Position',[0   130 120 20]);
      inp.probe.lengthWire.value          = uicontrol(uiPar{:},'String','','style','edit',       'Position',[120 130 70  20],'backgroundcolor','white','Callback',@(src,evt)obj.get_probe_length_wire);
      inp.probe.lengthWire.SIconversion   = 1e-2;
      inp.probe.radiusWire.text           = uicontrol(uiPar{:},'String','cyl/wire radius [cm]',  'Position',[0   110 120 20]);
      inp.probe.radiusWire.value          = uicontrol(uiPar{:},'String','',                      'Position',[120 110 70  20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_probe_radius_wire);
      inp.probe.radiusWire.SIconversion   = 1e-2;
      inp.probe.biasCurrent.text          = uicontrol(uiPar{:},'String','bias current [uA]',     'Position',[0   090 120 20]);
      inp.probe.biasCurrent.value         = uicontrol(uiPar{:},'String','0','style','edit',      'Position',[120 090 70  20],'backgroundcolor','white','Callback',@(src,evt)obj.get_probe_bias);
      inp.probe.biasCurrent.SIconversion  = 1e-6;
      %% initialize parameters menu
      inp.factorUvText                   = uicontrol('Parent',hp,'String','UV factor',             'Position',[0   70 60 20]);
      inp.factorUvValue                  = uicontrol('Parent',hp,'String',num2str(obj.InputParameters.factorUV),'Position',[70  70 100 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_factor_uv);
      inp.rSunText                       = uicontrol('Parent',hp,'String','Rsun [AU]',             'Position',[0   50 60 20]);
      inp.rSunValue                      = uicontrol('Parent',hp,'String',num2str(obj.InputParameters.rSunAU),              'Position',[70  50 100 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_distance_to_sun_au);
      inp.vectorUText                    = uicontrol('Parent',hp,'String','U [V]',                 'Position',[0   30 60 20]);
      inp.vectorUValue                   = uicontrol('Parent',hp,'String',num2str(obj.InputParameters.vectorUString),                      'Position',[70  30 100 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_u_interval);
      inp.update                         = uicontrol('Parent',hp,'String','Update',                'Position',[0   0 60 30],'Callback',@(src,evt)obj.calculate_ui);
      inp.reset                          = uicontrol('Parent',hp,'String','Reset',                 'Position',[70  0 60 30]); % TODO
      %% initialize s/c menu
      colPanelBg = [.95 1 1];
      hsc = uipanel('Title','Spacecraft','FontSize',12,'BackgroundColor',colPanelBg,'Position',[.7 .39 .3 .35]);
      popuptxt = obj.popup_list(obj.SpacecraftList);
      uiPar = {'Parent',hsc,'backgroundcolor',colPanelBg,'style','text'};
      inp.doModelSc                          = uicontrol(uiPar{:},'String','Model spacecraft','Value',0,      'Position',[0   205 120 25],'style','radio');
      inp.sc.name.text                       = uicontrol(uiPar{:},'String','spacecraft',                      'Position',[0   180 60  20]);
      inp.sc.name.value                      = uicontrol(uiPar{:},'String',popuptxt,                          'Position',[60  180 150 20],'style','popup','Callback',@(src,evt)obj.set_sc_model(src,evt));
      inp.sc.surface.text                    = uicontrol(uiPar{:},'String','surface',                         'Position',[0   160 60  20]);
      inp.sc.surface.value                   = uicontrol(uiPar{:},'String',[{'user defined'} lp.photocurrent],'Position',[60  160 150 20],'style','popup');
      inp.sc.areaTotal.text                  = uicontrol(uiPar{:},'String','Total area [m2]',                 'Position',[0   140 120 20]);
      inp.sc.areaTotal.value                 = uicontrol(uiPar{:},'String',num2str(0),'style','edit',         'Position',[120 140 50  20]);
      inp.sc.areaSunlit.text                 = uicontrol(uiPar{:},'String','Sunlit area [m2]',                'Position',[0   120 120 20],'Callback',@(src,evt)obj.get_sc_area_sunlit);
      inp.sc.areaSunlit.value                = uicontrol(uiPar{:},'String',num2str(0),'style','edit',         'Position',[120 120 50  20],'backgroundcolor','white');
      inp.sc.areaSunlitGuard.text            = uicontrol(uiPar{:},'String','Sunlit guard area [m2]',          'Position',[0   100 120  20],'Tooltipstring','Cross section area of pucks and guards, assuming similar photoelectron emission as antenna');
      inp.sc.areaSunlitGuard.value           = uicontrol(uiPar{:},'String','0','style','edit',                'Position',[120 100  50  20],'backgroundcolor','white');
      inp.sc.probeRefPotVsSatPot.text        = uicontrol(uiPar{:},'String','Probe refpot/scpot',              'Position',[0   80 120  20],'Tooltipstring','The ratio between the probe reference potential and satellite potential');
      inp.sc.probeRefPotVsSatPot.value       = uicontrol(uiPar{:},'String',num2str(0),                        'Position',[120 80  50  20],'style','edit','backgroundcolor','white');
      inp.sc.nProbes.text                    = uicontrol(uiPar{:},'String','Number of probes',                'Position',[0   60 120  20]);
      inp.sc.nProbes.value                   = uicontrol(uiPar{:},'String',num2str(0),                        'Position',[120 60  50  20],'style','edit','backgroundcolor','white');
      inp.sc.probeDistanceToSpacecraft.text  = uicontrol(uiPar{:},'String','distance probe-sc [m]',           'Position',[0   40 120  20]);
      inp.sc.probeDistanceToSpacecraft.value = uicontrol(uiPar{:},'String',num2str(0),                        'Position',[120 40  50  20],'style','edit','backgroundcolor','white');
      %% initialize plasma menu
      colPanelBg = [1 1 .95];
      hpl= uipanel('Title','Plasma','FontSize',12,'BackgroundColor',colPanelBg,'Position',[.7 .74 .3 .2]);
      popuptxt = obj.popup_list(obj.PlasmaList);
      uiPar = {'Parent',hpl,'backgroundcolor',colPanelBg,'style','text'};
      inp.plasma.typeText = uicontrol(uiPar{:},'String','model',      'Position',[0  100 60  20],'style','text');
      inp.plasma.typeValue= uicontrol(uiPar{:},'String',popuptxt,     'Position',[60 100 130 20],'style','popup','backgroundcolor','white','Callback',@(src,evt)obj.set_plasma_model(src,evt));
      inp.plasma.nString  = uicontrol(uiPar{:},'String','Ne [cc]',    'Position',[0 0 80 20]);
      inp.plasma.nValue   = uicontrol(uiPar{:},'String','',           'Position',[80 0 90 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_plasma_n);
      inp.plasma.TString  = uicontrol(uiPar{:},'String','T [eV]',     'Position',[0 20 80 20]);
      inp.plasma.TValue   = uicontrol(uiPar{:},'String','',           'Position',[80 20 90 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_plasma_t);
      inp.plasma.mpString = uicontrol(uiPar{:},'String','m [mp],0=me','Position',[0 40 80 20]);
      inp.plasma.mpValue  = uicontrol(uiPar{:},'String','',           'Position',[80 40 90 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_plasma_mp);
      inp.plasma.qeString = uicontrol(uiPar{:},'String','q [e]',      'Position',[0 60 80 20]);
      inp.plasma.qeValue  = uicontrol(uiPar{:},'String','',           'Position',[80 60 90 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_plasma_qe);
      inp.plasma.vString  = uicontrol(uiPar{:},'String','Vsc [km/s]', 'Position',[0 80 80 20]);
      inp.plasma.vValue   = uicontrol(uiPar{:},'String','',           'Position',[80 80 90 20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.get_plasma_v);
      %% initialize plot menu
      inp.toppanel.plotTypes = {'-','Resistance','Satellite IU','Antenna noise'};
      popuptxt = inp.toppanel.plotTypes;
      hpl               = uipanel('Title','Top panel','FontSize',12,'BackgroundColor',[1 1 .95],'Position',[.7 .94 .3 .06]);
      inp.toppanel.plot = uicontrol('Parent',hpl,'String',popuptxt,...
        'Position',[0 0 150 25],'style','popup','backgroundcolor','white');%TODO
      ud.inp            = inp;
      obj.figHandle     = figH;
      obj.UserData      = ud;

    end

    function set_plasma_model(obj,varargin)
      if nargin == 2 && any(strcmpi('user defined',varargin{1}))  %set_plasma_model(obj,'user defined')
        obj.plasmaUsed = 1;
        set(obj.UserData.inp.plasma.typeValue,'Value',1);
        return;
      elseif nargin == 2 && isnumeric(varargin{1}) %set_plasma_model(obj,numberInPlasmaList)
        idPlasma = varargin{1};
      elseif nargin == 3 %set_plasma_type(obj,hEvent,event)
        hEvent = varargin{1};
        event = varargin{2};
        disp(event);
        idPlasma = get(hEvent,'Value');
      else
        error('lp.ui.set_plasma_model unknown input');
      end
      set(obj.UserData.inp.plasma.typeValue,'Value',idPlasma);
      obj.plasmaUsed = idPlasma;
      plasma = obj.PlasmaList(idPlasma);
      set(obj.UserData.inp.plasma.nValue,'String',obj.field_to_vector_string(plasma,'n',1e-6));
      set(obj.UserData.inp.plasma.mpValue,'String',obj.field_to_vector_string(plasma,'mp'));
      set(obj.UserData.inp.plasma.TValue,'String',obj.field_to_vector_string(plasma,'T'));
      set(obj.UserData.inp.plasma.qeValue,'String',obj.field_to_vector_string(plasma,'qe'));
      set(obj.UserData.inp.plasma.vValue,'String',obj.field_to_vector_string(plasma,'v',1e-3));
    end
    function get_plasma_qe(obj)
      qeStr = get(obj.UserData.inp.plasma.qeValue,'String'); % in [e]
      if isempty(qeStr), qeStr = '-1'; end
      qe = eval(['[' qeStr ']']);
      obj.set_user_defined_if_plasma_changes('qe',qe)
    end
    function get_plasma_mp(obj)
      mpStr = get(obj.UserData.inp.plasma.mpValue,'String'); % in [e]
      if isempty(mpStr), mpStr = '-1'; end
      mp = eval(['[' mpStr ']']);
      obj.set_user_defined_if_plasma_changes('mp',mp)
    end
    function get_plasma_n(obj)
      nStr = get(obj.UserData.inp.plasma.nValue,'String'); % in [e]
      if isempty(nStr), nStr = '-1'; end
      n = eval(['[' nStr ']'])*1e6;
      obj.set_user_defined_if_plasma_changes('n',n)
    end
    function get_plasma_t(obj)
      tStr = get(obj.UserData.inp.plasma.TValue,'String'); % in [e]
      if isempty(tStr), tStr = '-1'; end
      T = eval(['[' tStr ']']);
      obj.set_user_defined_if_plasma_changes('T',T)
    end
    function get_plasma_v(obj)
      vStr = get(obj.UserData.inp.plasma.vValue,'String'); % in [e]
      if isempty(vStr), vStr = '-1'; end
      v = eval(['[' vStr ']'])*1e3; % [km/s] > [m/s]
      obj.set_user_defined_if_plasma_changes('v',v)
    end

    function set_sc_model(obj,varargin)
      if nargin == 2 && any(strcmpi('user defined',varargin{1}))  %set_sc_model(obj,'user defined')
        obj.spacecraftUsed = 1;
        set(obj.UserData.inp.spacecraft.typeValue,'Value',1);
        return;
      elseif nargin == 2 && isnumeric(varargin{1}) %set_sc_model(obj,numberInSpacecraftList)
        idSpacecraft = varargin{1};
      elseif nargin == 3 %set_plasma_type(obj,hEvent,event)
        hEvent = varargin{1};
        event = varargin{2};
        disp(event);
        idSpacecraft = get(hEvent,'Value');
      else
        error('lp.ui.set_sc_model unknown input');
      end
      set(obj.UserData.inp.sc.name.value,'Value',idSpacecraft);
      obj.spacecraftUsed = idSpacecraft;
      spacecraft = obj.SpacecraftList(idSpacecraft);
      fieldsToUpdate = {'areaTotal','areaSunlit','areaSunlitGuard','probeRefPotVsSatPot','nProbes','probeDistanceToSpacecraft'};
      obj.UserData.inp.sc = obj.update_input_fields(obj.UserData.inp.sc,fieldsToUpdate,spacecraft);
    end
    function get_sc_area_sunlit(obj)
      areaSunlitStr = get(obj.UserData.inp.sc.areaSunlitValue,'String'); % in m^2
      if isempty(areaSunlitStr), areaSunlitStr = '0'; end
      areaSunlit = str2double(areaSunlitStr)*1e-2;
      obj.set_user_defined_if_sc_changes('areaSunlit',areaSunlit)
    end
    function get_spacecraft_area_total(obj)
      areaTotal = get(obj.UserData.inp.sc.areaTotal.value,'String'); % in cm
      if isempty(areaTotal), areaTotal = 0; end
      lengthWire = str2double(areaTotal)*1e-2;
      obj.set_user_defined_if_probe_changes('lengthWire',lengthWire)
    end
    function get_factor_uv(obj)
      factorUvString = get(obj.UserData.inp.factorUvValue,'String'); % in cm
      if isempty(factorUvString)
        factorUvString = '1';
        set(obj.UserData.inp.factorUvValue,'String','1')
      end
      obj.InputParameters.factorUV = str2double(factorUvString);
    end
    function get_distance_to_sun_au(obj)
      rSunString = get(obj.UserData.inp.rSunValue,'String'); % in cm
      if isempty(rSunString)
        rSunString = '1';
        set(obj.UserData.inp.rSunValue,'String','1')
      end
      obj.InputParameters.rSunAU = str2double(rSunString);
    end
    function get_u_interval(obj)
      vectorUString = get(obj.UserData.inp.vectorUValue,'String'); % in cm
      obj.InputParameters.vectorUonlyMinMaxGiven = false;
      if isempty(vectorUString)
        vectorUString = '-5:10';
        set(obj.UserData.inp.vectorUValue,'String',vectorUString);
      end
      obj.InputParameters.vectorUString = vectorUString;
      obj.InputParameters.vectorU = eval(['[' vectorUString ']']);
      if numel(obj.InputParameters.vectorU) == 2 % only min and max are given
        obj.InputParameters.vectorU = obj.InputParameters.vectorU(1)*[1 0.75 0.5 0.25 0] +...
          obj.InputParameters.vectorU(2)*[0 0.25 0.5 0.75 1];
        obj.InputParameters.vectorUonlyMinMaxGiven = true;
      end
    end

    function set_probe_model(obj,varargin)
      if nargin == 2 && any(strcmpi('user defined',varargin{1}))  %set_probe_model(obj,'user defined')
        idProbe = 1;
      elseif nargin == 2 && isnumeric(varargin{1}) %set_probe_model(obj,numberOfProbe)
        idProbe = varargin{1};
      elseif nargin == 3 %set_probe_model(obj,hEvent,event)
        hEvent = varargin{1};
        idProbe = get(hEvent,'Value');
      end
      obj.probeUsed = idProbe;
      obj.update_probe_area_total_vs_sunlit;
      obj.set_probe_surface;
      obj.set_probe_radius_sphere(obj.ProbeList(idProbe).radiusSphere);
      obj.set_probe_radius_wire(  obj.ProbeList(idProbe).radiusWire);
      obj.set_probe_length_wire(  obj.ProbeList(idProbe).lengthWire);
      set(obj.UserData.inp.probe.type.value,'Value',idProbe);
    end
    function set_probe_surface(obj,varargin)
      if nargin == 1
        probeParameters = obj.ProbeList(obj.probeUsed);
        idProbeSurface = find(strcmp(probeParameters.surface,obj.ProbeSurfaceList));
        if idProbeSurface
          obj.set_probe_surface(idProbeSurface);
        else
          irf.log('critical',[surface '''' probeParameters.surface ''' is unknown by lp.photocurrent.']);
        end
      elseif nargin == 2 && isnumeric(varargin{1}) %set_probe_model(obj,numberOfSurface)
        obj.probeSurfaceUsed = varargin{1};
        set(obj.UserData.inp.probe.surface.value,'value',obj.probeSurfaceUsed);
        set(obj.UserData.inp.probe.surfacePhotoemission.value,...
          'String',num2str(obj.ProbeSurfacePhotoemission(obj.probeSurfaceUsed)*1e6));
      elseif nargin == 3 %set_probe_model(obj,hEvent,event)
        hEvent = varargin{1};
        idProbeSurface = get(hEvent,'Value');
        obj.set_probe_surface(idProbeSurface);
      end
    end
    function get_probe_surface(obj)
      surfacePhotoemissionMicroA = obj.UserData.inp.probe.surfacePhotoemission.value.String; % in cm
      if isempty(surfacePhotoemissionMicroA), surfacePhotoemissionMicroA = '0'; end
      surfacePhotoemission = str2double(surfacePhotoemissionMicroA)*1e-6;
      if surfacePhotoemission ~= obj.ProbeSurfacePhotoemission(obj.probeSurfaceUsed)
        obj.ProbeSurfacePhotoemission(1) = surfacePhotoemission;
        obj.set_probe_surface(1);
      end
    end
    function set_probe_radius_sphere(obj,radiusSphere)
      set(obj.UserData.inp.probe.radiusSphere.value,...
        'String',num2str( radiusSphere/obj.UserData.inp.probe.radiusSphere.SIconversion,2 ));
      obj.update_probe_area_total_vs_sunlit;
    end
    function get_probe_radius_sphere(obj)
      radiusSphereCm = get(obj.UserData.inp.probe.radiusSphere.value,'String'); % in cm
      if isempty(radiusSphereCm), radiusSphereCm = 0; end
      radiusSphere = str2double(radiusSphereCm)*1e-2;
      obj.set_user_defined_if_probe_changes('radiusSphere',radiusSphere)
    end
    function get_probe_radius_wire(obj)
      radiusWireCm = get(obj.UserData.inp.probe.radiusWire.value,'String'); % in cm
      if isempty(radiusWireCm), radiusWireCm = 0; end
      radiusWire = eval(['[' radiusWireCm ']'])*1e-2;
      obj.set_user_defined_if_probe_changes('radiusWire',radiusWire)
    end
    function get_probe_length_wire(obj)
      lengthWireCm = get(obj.UserData.inp.probe.lengthWire.value,'String'); % in cm
      if isempty(lengthWireCm), lengthWireCm = 0; end
      lengthWire = str2double(lengthWireCm)*1e-2;
      obj.set_user_defined_if_probe_changes('lengthWire',lengthWire)
    end
    function set_probe_radius_wire(obj,radiusWire)
      set(obj.UserData.inp.probe.radiusWire.value,...
        'String',num2str(radiusWire/obj.UserData.inp.probe.radiusWire.SIconversion,2));
      obj.update_probe_area_total_vs_sunlit;
    end
    function set_probe_length_wire(obj,lengthWire)
      set(obj.UserData.inp.probe.lengthWire.value,...
        'String',num2str(lengthWire/obj.UserData.inp.probe.lengthWire.SIconversion,3));
      obj.update_probe_area_total_vs_sunlit;
    end
    function get_probe_bias(obj) % input in [A]
      biasCurrenMicroA = get(obj.UserData.inp.probe.biasCurrent.value,'String'); % in [uA]
      obj.InputParameters.biasCurrent = str2double(biasCurrenMicroA) / 1e6;
    end
    function update_probe_area_total_vs_sunlit(obj)
      set(obj.UserData.inp.probe.areaTotalVsSunlit.value,...
        'String',num2str(obj.ProbeList(obj.probeUsed).Area.totalVsSunlit,3));
    end

    function calculate_ui(obj)
      doModelSpacecraft = obj.UserData.inp.doModelSc.Value;
      if doModelSpacecraft
        obj.get_probe_and_spacecraft_current;
      else
        obj.get_probe_current;
      end
      obj.plot_ui;
    end
    function get_probe_current(obj)
      ProbeToUse = obj.ProbeList(obj.probeUsed);
      ProbeToUse.surfacePhotoemission = obj.ProbeSurfacePhotoemission(obj.probeSurfaceUsed);
      slopeChangeMax=1.3; % how big maximum relative slope change is allowed
      slopeChangeMin=1+(slopeChangeMax-1)/4;
      vecU = obj.InputParameters.vectorU;
      doIterationToFindUI = true;
      while doIterationToFindUI
        jProbe = lp.current(ProbeToUse,vecU,...
          obj.InputParameters.rSunAU,...
          obj.InputParameters.factorUV,...
          obj.PlasmaList(obj.plasmaUsed));
        I = jProbe.total;
        [U__dUdI,dU,dUdI] = obj.dUdI_(vecU,I);
        if obj.InputParameters.vectorUonlyMinMaxGiven
          doIterationToRemoveMins = true;
          vecUnew = vecU;
          while doIterationToRemoveMins
            [U__dUdI,dU,dUdI] = obj.dUdI_(vecUnew,I);
            slopeChanges = abs(dUdI(2:end)./dUdI(1:end-1));
            ii = find(slopeChanges<1);
            slopeChanges(ii) = 1./slopeChanges(ii);
            iMin = find(slopeChanges < slopeChangeMin)+1; % find points where slope changes little in both directions from the point
            if any(iMin)
              vecUnew(iMin(1:2:end)) = [];
              I(iMin(1:2:end)) =[];
            else
              doIterationToRemoveMins = false;
            end
          end
          iMax = find((slopeChanges > slopeChangeMax) & ((dU(1:end-1) > 0.2) & (dU(2:end) > 0.2)));
          if any(iMax)
            for j=numel(iMax):-1:1
              vecUnew = [vecUnew(1:iMax(j)) ...
                [vecUnew(iMax(j))   vecUnew(iMax(j)+1)]*[0.7 0.3;0.3 0.7]...
                [vecUnew(iMax(j)+1) vecUnew(iMax(j)+2)]*[0.7 0.3;0.3 0.7]...
                vecUnew(iMax(j)+2:end)];
            end
          end
          if isempty(iMax)
            doIterationToFindUI = false;
          else
            vecU = vecUnew;
          end
        else
          doIterationToFindUI = false;
        end
      end
      obj.Output.J = jProbe;
      obj.Output.U = vecU;
      obj.Output.U__dUdI = U__dUdI;
      obj.Output.dUdI = dUdI;
    end
    function get_probe_and_spacecraft_current(obj)
      obj.get_probe_current;
      Upot   = obj.Output.U;
      jProbe = obj.Output.J.total;
      obj.Output.J.sc=lp.current(obj.SpacecraftList(obj.spacecraftUsed),...
        Upot,...
        obj.InputParameters.rSunAU,...
        obj.InputParameters.factorUV,...
        obj.PlasmaList(obj.plasmaUsed));
      obj.Output.UI.dUdIsat=gradient(Upot,obj.Output.J.sc.total);
      % reduce probe curve to reasonable number of points (derivative does
      % not change more than 10% between points
      Ibias = jProbe; % bias current to probe measured with respect to sc
      Ubias = Upot;

      %
      Probe = obj.ProbeList(obj.probeUsed);
      Sc    = obj.SpacecraftList(obj.spacecraftUsed);
      antena_guard_area_factor=(Probe.Area.sunlit+Sc.areaSunlitGuard/Sc.nProbes)/Probe.Area.sunlit;
      Iprobe=Ibias*antena_guard_area_factor;
      % zero approximation
      Isat=Sc.nProbes*Iprobe;
      Usatsweep=interp1(obj.Output.J.sc.total,Upot,Isat/5,'nearest'); % floating potential of sc during sweep, assuming that only 1/5 of bias current electrons esacpe to space
      Isat_probe_photoelectrons=Isat*0; % electrons from probes hitting spacecraft
      Uprobe2plasma=zeros(size(Ibias)); % initialize
      Uproberefsweep=zeros(size(Ibias)); % initialize
      Jprobephotoreturn=zeros(size(Ibias)); % initialize
      FitError=zeros(size(Ibias)); % initialize
      refpotvec=Ubias(1:4:end)*Sc.probeRefPotVsSatPot;
      [probepotgrid,refpotgrid] = meshgrid(Ubias,refpotvec);
      Jprobephotogrid=lp.current(Probe,probepotgrid-refpotgrid,obj.InputParameters.rSunAU,obj.InputParameters.factorUV,[]);
      Jprobephotoescapingscgrid=lp.current(Probe,probepotgrid,obj.InputParameters.rSunAU,obj.InputParameters.factorUV,[]);
      Jprobeplasmagrid=lp.current(Probe,probepotgrid,obj.InputParameters.rSunAU,0,obj.PlasmaList(obj.plasmaUsed));
      jProbegrid=Jprobeplasmagrid.total+Jprobephotogrid.total;
      jProbephotoe2scgrid=Jprobephotogrid.total-Jprobephotoescapingscgrid.total;
      for ii=1:numel(Iprobe)
        % plasma current with UV factor zero
        satpot=Usatsweep(ii);
        uprobe=interp1(jProbe,Upot,Ibias(ii))+satpot*Sc.probeRefPotVsSatPot;
        ibias=Ibias(ii);
        jj=1;
        err=1;
        while err > .001 && jj<10
          refpot=satpot*Sc.probeRefPotVsSatPot;
          if isnan(refpot), break; end
          uprobe=interp1(interp1(refpotvec,jProbegrid,refpot,'linear','extrap'),Ubias,ibias);
          if isnan(uprobe); uprobe=interp1(jProbe,Upot,Ibias(ii))+satpot*Sc.probeRefPotVsSatPot;end
          J_probe_photoe2sc=interp2(probepotgrid,refpotgrid,jProbephotoe2scgrid,uprobe,refpot);
          new_Isat_probe_photoelectrons=-J_probe_photoe2sc*antena_guard_area_factor*Sc.nProbes;
          err=abs(new_Isat_probe_photoelectrons-Isat_probe_photoelectrons(ii))/abs(new_Isat_probe_photoelectrons);
          if isnan(err), err=0; end
          Isat_probe_photoelectrons(ii)=new_Isat_probe_photoelectrons;
          satpot=interp1(obj.Output.J.sc.total,Upot,-Isat(ii)-Isat_probe_photoelectrons(ii));
          jj=jj+1;
        end
        FitError(ii)=err;
        Jprobephotoreturn(ii)=J_probe_photoe2sc;
        Usatsweep(ii)=satpot;
        Uprobe2plasma(ii)=uprobe;
        Uproberefsweep(ii)=refpot;
      end

      obj.Output.UI.Jprobephotoreturn= Jprobephotoreturn;
      obj.Output.UI.Ubias            = Ubias;
      obj.Output.UI.Usatsweep        = Usatsweep;
      obj.Output.UI.Ibias            = Ibias;
      obj.Output.UI.Uprobe2plasma    = Uprobe2plasma;
      obj.Output.UI.Uprobe2sc        = Uprobe2plasma-Usatsweep;
      obj.Output.UI.Uprobe2refpot    = Uprobe2plasma-Uproberefsweep;
      obj.Output.UI.dUdI_probe2plasma= gradient(Uprobe2plasma,Ibias);
      obj.Output.UI.dUdI_probe2sc    = gradient(obj.Output.UI.Uprobe2sc,Ibias);
      obj.Output.UI.dUdI_probe2refpot= gradient(obj.Output.UI.Uprobe2refpot,Ibias);
      obj.Output.UI.dUdI             = obj.Output.UI.dUdI_probe2refpot;
      obj.Output.UI.Upot             = obj.Output.UI.Uprobe2refpot;
      obj.Output.UI.J_probe          = Ibias;
    end
    function plot_ui(obj)
      % Bottom axes
      obj.plot_UI_probe;

      %Top axes
      toppanelPlotType = get(obj.UserData.inp.toppanel.plot,'Value');
      switch toppanelPlotType
        case 1 % Nothing
          cla(obj.Axes.top);
        case 2 % Probe Resistance
          obj.plot_resistance;
          linkaxes([obj.Axes.bottom, obj.Axes.top],'x');
      end

      % Information panel
      obj.plot_UI_information_panel;
    end
    function plot_resistance(obj)
      h=obj.Axes.top;
      U = obj.Output.U;
      I = obj.Output.J.total;
      vecU = (U(1:end-1) + U(2:end) ) / 2;
      dU = U(2:end) - U(1:end-1);
      dI = I(2:end) - I(1:end-1);
      dUdI = dU./dI;
      plot(h,vecU,dUdI,'k');
      set(h,'xlim',[min(vecU) max(vecU)]);
      grid(h,'on');
      xlabel(h,'U [V]');
      ylabel(h,'dU/dI [\Omega]');
      if 0%ud.flag_use_sc, % add probe resistance wrt plasma and s/c
        hold(h(2),'on'); %#ok<UNRCH>
        plot(h(2),Uprobe2plasma,dUdI_probe2plasma,'r','linewidth',1.5);
        plot(h(2),Uprobe2sc,dUdI_probe2sc,'b','linewidth',1.5);
        hold(h(2),'off');
      end
      axis(h,'auto y');
      set(h,'yscale','log')
    end
    function plot_UI_probe(obj)
      % plot IU curve
      %info_txt='';
      h=obj.Axes.bottom;
      doModelSpacecraft = obj.UserData.inp.doModelSc.Value;
      vecU = obj.Output.U;
      J = obj.Output.J;
      biasCurrentA = obj.InputParameters.biasCurrent;
      biasCurrentMicroA = 1e6*biasCurrentA;
      flagBias = biasCurrentMicroA ~= 0 && biasCurrentA>min(J.total) && -biasCurrentA<max(J.total);

      plot(h,vecU, J.total*1e6,'k');
      set(h,'xlim',[min(vecU) max(vecU)]);
      grid(h,'on');
      xlabel(h,'U [V]');
      ylabel(h,'I [\mu A]');

      % Add photoelectron current
      hold(h,'on');
      if doModelSpacecraft
        plot(h,obj.Output.UI.Ubias,obj.Output.UI.Ibias*1e6,'k','linewidth',0.5);
        plot(h,obj.Output.UI.Uprobe2plasma,obj.Output.UI.Ibias*1e6,'r.','linewidth',1.5);
        plot(h,obj.Output.UI.Uprobe2sc,obj.Output.UI.Ibias*1e6,'b','linewidth',1.5);
        plot(h,obj.Output.UI.Uprobe2sc,obj.Output.UI.Jprobephotoreturn*1e6,'color',[0 0.5 0],'linewidth',1.5);
        plot(h,obj.Output.UI.Usatsweep,obj.Output.UI.Ibias*1e6,'color',[0.5 0 0.5],'linewidth',1);
        irf_legend(h,'probe to plasma',[1 1.01],'color','r');
        irf_legend(h,'probe to reference',[1 1.05],'color','k');
        irf_legend(h,'probe to spacecraft (bias)',[0.02 1.01],'color','b');
        irf_legend(h,'probe photo e- to s/c',[0.02 1.05],'color',[0 0.5 0]);
        irf_legend(h,'Satellite potential',[0.98 0.03],'color',[0.5 0 0.5]);
        axis(h,'auto x');
        if flagBias % draw bias current line
          plot(h,[obj.Output.UI.Uprobe2sc(1) obj.Output.UI.Uprobe2plasma(end)],biasCurrentMicroA.*[-1 -1],'k-.','linewidth',0.5);
          text(obj.Output.UI.Uprobe2sc(1),-biasCurrentMicroA,'bias','parent',h,'horizontalalignment','left','verticalalignment','bottom');
          %flag_add_bias_point_values=1; % TODO
        end
      else
        plot(h,vecU,J.photo*1e6,'r','linewidth',0.5);
        irf_legend(h,'    total',      [0.98 0.03],'color','k');
        irf_legend(h,' photoelectrons',[0.98 0.08],'color','r');
        clr=[0.5 0 0; 0 0.5 0; 0 0 0.5];
        for ii=1:length(J.plasma)
          plot(h,vecU,J.plasma{ii}*1e6,'linewidth',.5,'color',clr(:,ii));
          irf_legend(h,['plasma ' num2str(ii)],[0.98 0.08+ii*0.05],'color',clr(:,ii));
        end
        if flagBias % draw bias current line
          plot(h,[vecU(1) vecU(end)],biasCurrentMicroA.*[-1 -1],'k-.','linewidth',0.5);
          text(vecU(1),-biasCurrentMicroA,'-bias','parent',h,'horizontalalignment','left','verticalalignment','bottom');
        end
      end
      hold(h,'off');

    end
    function plot_UI_information_panel(obj)
      % prepare informationa panel

      % Calculate all required values
      dUdI = obj.Output.dUdI;
      U__dUdI = obj.Output.U__dUdI;
      vecU = obj.Output.U;
      J    = obj.Output.J;
      biasCurrentA = obj.InputParameters.biasCurrent;
      biasCurrentMicroA = 1e6*biasCurrentA;
      flagBias = biasCurrentMicroA ~= 0 && biasCurrentA>min(J.total) && -biasCurrentA<max(J.total);
      doModelSpacecraft = obj.UserData.inp.doModelSc.Value;
      probe = obj.ProbeList(obj.probeUsed);

      Rmin = min(abs(dUdI)); % minimum resistance
      fcr=1/2/pi/Rmin/probe.capacitance;
      disp(['Rmin=' num2str(Rmin,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, f_{CR}=' num2str(fcr,3) 'Hz.']);
      InfoTxt = struct();
      if doModelSpacecraft
        InfoTxt.probeToPlasma =['probe to plasma Rmin=' num2str(min(abs(obj.Output.UI.dUdI_probe2plasma)),3) ' Ohm'];
        InfoTxt.probeToRef    =['probe to reference Rmin=' num2str(min(abs(obj.Output.UI.dUdI)),3) ' Ohm'];
        InfoTxt.probeToSc     =['probe to spacecraft Rmin=' num2str(min(abs(obj.Output.UI.dUdI_probe2sc)),3) ' Ohm'];
      else
        InfoTxt.probe=['Probe at minimum R: R =' num2str(Rmin,3) ' Ohm,'...
          ' C =' num2str(probe.capacitance*1e12,3) 'pF,' ...
          ' fcr =' num2str(fcr,3) 'Hz.'];
      end

      if min(J.total)<0 && max(J.total)>0 % display information on Ufloat
        Ufloat=interp1(J.total,vecU,0);    % floating potential
        ii=isfinite(U__dUdI);
        Rfloat=interp1(U__dUdI(ii),dUdI(ii),Ufloat);
        fcr=1/2/pi/Rfloat/probe.capacitance;
        InfoTxt.probeFloat =['Floating probe: U =' num2str(Ufloat,3) 'V,' ...
          ' R = ' num2str(Rfloat,3) ' Ohm,' ...
          ' C =' num2str(probe.capacitance*1e12,3) 'pF,' ...
          ' fcr=' num2str(fcr,3) 'Hz.'];
        disp(['Probe: Ufloat=' num2str(Ufloat,3) ' V, Rfloat=' num2str(Rfloat,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.']);
      end
      if flagBias %flag_add_bias_point_values,
        biasCurrentA = obj.InputParameters.biasCurrent;
        Ubias=interp1(J.total,vecU,-biasCurrentA); % floating potential
        ii=isfinite(U__dUdI);
        Rbias=interp1(U__dUdI(ii),dUdI(ii),Ubias);
        fcr=1/2/pi/Rbias/probe.capacitance;
        InfoTxt.probeBias =['Biased probe: U =' num2str(Ubias,3) 'V,' ...
          ' R = ' num2str(Rbias,3) ' Ohm,' ...
          ' C =' num2str(probe.capacitance*1e12,3) 'pF,' ...
          ' fcr=' num2str(fcr,3) 'Hz.'];
        disp(['Rbias=' num2str(Rbias,3) ' Ohm, C=' num2str(probe.capacitance*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.']);
        if doModelSpacecraft
          Uscbias=interp1(J.total,obj.Output.UI.Usatsweep,-biasCurrentA); % floating potential
          ii=isfinite(vecU);
          Rscbias=interp1(vecU(ii),obj.Output.UI.dUdIsat(ii),Uscbias);
          cSat = irf_estimate('capacitance_sphere',sqrt(obj.SpacecraftList(obj.spacecraftUsed).areaTotal/(4*pi)));
          fcr=1/2/pi/Rscbias/cSat;
          disp(['Spacecraft (biased): Rbias=' num2str(Rbias,3) ' Ohm, C=' num2str(cSat*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.']);
          InfoTxt.biasedSpacecraft =['Spacecraft (biased): Ubias=' num2str(Uscbias,3)  ', Rbias=' num2str(Rscbias,3) 'Ohm, C=' num2str(cSat*1e12,3) 'pF, fcr=' num2str(fcr,3) 'Hz.'];
          Usp=-interp1(obj.Output.UI.Ibias,obj.Output.UI.Uprobe2sc,-biasCurrentA); % floating potential
          disp(['Spacecraft to Probe: Usp=' num2str(Usp,3) 'V.']);
          InfoTxt.spacecraftToProbe = ['Spacecraft to Probe: Usp=' num2str(Usp,3)  ' V.'];
        end
      end
      if 0%ud.flag_use_sc && min(ud.I_sc)<0 && max(ud.I_sc)>0, % display information on Ufloat
        Uscfloat=interp1(ud.I_sc,ud.U_sc,0);  %#ok<UNRCH> % floating potential
        ii=isfinite(ud.U_sc);Rscfloat=interp1(ud.U_sc(ii),ud.dUdI_sc(ii),Uscfloat);
        info_txt=[info_txt '\newline Spacecraft: Ufloat=' num2str(Uscfloat,3) 'V, Rfloat= ' num2str(Rscfloat,3) ' Ohm, C=' num2str(ud.sc.capacitance*1e12,3) 'pF'];
        disp(['Spacecraft: Ufloat=' num2str(Uscfloat,3) ' V, Rfloat=' num2str(Rscfloat,3) ' Ohm, C=' num2str(ud.sc.capacitance*1e12,3) 'pF']);
      end
      if 0%ud.UV_factor>0,                                     % display photoelectron saturation current
        infoTxt=[info_txt '\newline Probe: photo e- Io = ' num2str(ud.UV_factor*lp.photocurrent(1,-1,ud.R_sun,ud.probe.surface)*1e6,3) '[\mu A/m^2]'];  %#ok<UNRCH>
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

      % Write out all required output
      fieldsInfo = fieldnames(InfoTxt);
      infoText = '';
      for iField = 1:numel(fieldsInfo)
        infoText = [infoText InfoTxt.(fieldsInfo{iField}) '\newline']; %#ok<AGROW>
      end
      set(obj.Axes.infoTextHandle,'string',infoText,'fontsize',11);
    end
    function set_user_defined_if_plasma_changes(obj,idString,vector)
      if (obj.plasmaUsed ~= 1) ...
          && numel(vector) == numel(obj.PlasmaList(obj.plasmaUsed).(idString)) ...
          && all(vector == obj.PlasmaList(obj.plasmaUsed).(idString))
        return
      else
        obj.PlasmaList(1) = obj.PlasmaList(obj.plasmaUsed);
        obj.PlasmaList(1).(idString) = vector;
        obj.set_plasma_model('user defined');
      end
    end
    function set_user_defined_if_probe_changes(obj,idString,vector)
      if (obj.probeUsed ~= 1) ...
          && numel(vector) == numel(obj.ProbeList(obj.probeUsed).(idString)) ...
          && all(vector == obj.ProbeList(obj.probeUsed).(idString))
        return
      else
        obj.ProbeList(1) = obj.ProbeList(obj.probeUsed);
        obj.ProbeList(1).(idString) = vector;
        obj.update_probe_area_total_vs_sunlit;
        obj.set_probe_model('user defined');
      end
    end
    function set_user_defined_if_sc_changes(obj,idString,vector)
      if (obj.spacecraftUsed ~= 1) ...
          && numel(vector) == numel(obj.SpacecraftList(obj.spacecraftUsed).(idString)) ...
          && all(vector == obj.SpacecraftList(obj.spacecraftUsed).(idString))
        return
      else
        obj.SpacecraftList(1) = obj.SpacecraftList(obj.spacecraftUsed);
        obj.SpacecratList(1).(idString) = vector;
        obj.update_sc_area_total_vs_sunlit;
        obj.set_sc_model('user defined');
      end
    end
  end
  methods (Static)
    function popupText = popup_list(inp)
      if iscell(inp) && (numel(inp) > 0)
        popupText =inp{1};
        for ii=2:numel(inp)
          popupText(end+1:end+1+numel(inp{ii}))=['|' inp{ii}];
        end
      elseif numel(inp) > 0 && all(isprop(inp,'name'))
        popupText ='';
        for ii=1:numel(inp)
          probeName = inp(ii).name;
          if ~isempty(probeName)
            popupText(end+1:end+1+numel(probeName))=[probeName '|'];
          end
        end
        popupText(end) =[];
      end
    end
    function str = field_to_vector_string(o,field,multFactor)
      str = '';
      if nargin == 2
        multFactor = 1;
      end
      for ii = 1:numel(o)
        str = [str num2str(o(ii).(field)*multFactor) ' ']; %#ok<AGROW>
      end
      str(end) = [];
    end
    function OutputGui = update_input_fields(InputGui,fieldsToUpdate,InputObject)
      for iField = 1:numel(fieldsToUpdate)
        field = fieldsToUpdate{iField};
        if isfield(InputGui.(field),'SIconversion')
          SIconversion = InputGui.(field).('SIconversion');
        else
          SIconversion = 1;
        end
        set(InputGui.(field).value,'String',lp.ui.field_to_vector_string(InputObject,field,SIconversion));
      end
      OutputGui = InputGui;
    end
    function [U__dUdI,dU,dUdI] = dUdI_(U,I)
      U__dUdI = (U(1:end-1) + U(2:end) ) / 2;
      dU = U(2:end) - U(1:end-1);
      dI = I(2:end) - I(1:end-1);
      dUdI = dU./dI;
    end
  end
end

