classdef ui < handle
	%LP.UI User interface to +lp package
	%   LP.UI(spacecraft,probe,plasma,parameters)
	properties (SetAccess = protected)
		SpacecraftList % Nr 1 is always user defined
		spacecraftUsed
		ProbeList % Nr 1 is always user defined  
		probeUsed % 1..N from ProbeList
		PlasmaList
		plasmaUsed
		figHandle
		UserData % user data
		InputParameters
	end
	methods
		function obj = ui(varargin)
			%% Message shown once during the session
			message='You can anytime access all the results from get(gcf,''userdata'').';
			disp(message);
			
			%% Check input
			args=varargin;
			while ~isempty(args)
				if isa(args{1},'lp.spacecraft'),
					obj.SpacecraftList = args{1};
					args(1) =[];
				elseif isa(args{1},'lp.lprobe'),
					obj.ProbeList = args{1};
					args(1) =[];
				elseif isa(args{1},'lp.plasma')
					obj.PlasmaList = args{1};
					args(1) =[];
				else
					error('lp.ui input unknown!');
				end
			end
			if ~isa(obj.SpacecraftList,'lp.spacecraft'),
				obj.SpacecraftList = lp.default_spacecraft;
			end
			if ~isa(obj.ProbeList,'lp.probe'),
				obj.ProbeList = lp.default_lprobe;
			end
			if ~isa(obj.PlasmaList,'lp.plasma'),
				obj.PlasmaList = lp.default_plasma;
			end
			obj.SpacecraftList = [obj.SpacecraftList(1) obj.SpacecraftList];
			obj.SpacecraftList(1).name = 'user defined';
			obj.ProbeList = [obj.ProbeList(1) obj.ProbeList];
			obj.ProbeList(1).name = 'user defined';
			obj.PlasmaList = [obj.PlasmaList(1) obj.PlasmaList];
			obj.PlasmaList(1).name = 'user defined';
			%% Initialize IDE
			obj.new_ide();
			obj.spacecraftUsed = 2;
			obj.probeUsed = 2;
			obj.plasmaUsed = 2;
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
			ud.h(1)=axes('position',[0.1 0.3 0.5 0.3]); % [x y dx dy]
			ud.h(2)=axes('position',[0.1 0.67 0.5 0.3]); % [x y dx dy]
			ud.h(3)=axes('position',[0.1 0.0 0.5 0.13]); % [x y dx dy]
			linkaxes(ud.h(1:2),'x');
			axis(ud.h(3),'off');
			ud.ht=text(0,1,'','parent',ud.h(3));
			set(figH,'userdata',ud);
			%% initialize probe menu
			colPanelBg = [1 0.95 1];
			hp = uipanel('Title','Probe','FontSize',12,'BackgroundColor',colPanelBg,'Position',[.7 .0 .3 .39]);
			popuptxt = obj.popup_list(obj.ProbeList);
			inp.probe.typeText                   = uicontrol('Parent',hp,'String','type','style','text',   'Position',[0   230 60   20],'backgroundcolor',colPanelBg);
			inp.probe.typeValue                  = uicontrol('Parent',hp,'String',popuptxt,'style','popup','Position',[60  230 130  20],'backgroundcolor','white','Callback',@(src,evt)obj.set_probe_type(src,evt));
			inp.probe.srufaceText                = uicontrol('Parent',hp,'String','surface','style','text','Position',[0   210 60   20],'backgroundcolor',colPanelBg);
			inp.probe.surfaceValue               = uicontrol('Parent',hp,'String',[{'User defined'} lp.photocurrent],'Position',[60  210 130  20],'style','popup','backgroundcolor','white');
			inp.probe.total_vs_sunlit_area_text  = uicontrol('Parent',hp,'String','total/sunlit area',     'Position',[0   190 120 20],'style','text');
			inp.probe.total_vs_sunlit_area_value = uicontrol('Parent',hp,'String','','style','text',       'Position',[120 190 70  20],'backgroundcolor','white');
			inp.probe.radiusSphereText           = uicontrol('Parent',hp,'String','sphere radius [cm]',    'Position',[0   170 120 20]);
			inp.probe.radiusSphereValue          = uicontrol('Parent',hp,'String','','style','edit',       'Position',[120 170 70  20],'backgroundcolor','white','Callback',@(src,evt)obj.get_probe_radius_sphere);
			inp.probe.lengthWireText             = uicontrol('Parent',hp,'String','cyl/wire length [cm]',  'Position',[0   150 120 20]);
			inp.probe.lengthWireValue            = uicontrol('Parent',hp,'String','','style','edit',       'Position',[120 150 70  20],'backgroundcolor','white','Callback',@(src,evt)obj.set_probe_length_wire(src,evt));
			inp.probe.radiusWireText             = uicontrol('Parent',hp,'String','cyl/wire radius [cm]',  'Position',[0   130 120 20]);
			inp.probe.radiusWireValue            = uicontrol('Parent',hp,'String','',                      'Position',[120 130 70  20],'style','edit','backgroundcolor','white','Callback',@(src,evt)obj.set_probe_radius_wire(src,evt));
			inp.probe.biasCurrentText            = uicontrol('Parent',hp,'String','bias current [uA]',     'Position',[0   110 120 20]);
			inp.probe.biasCurrentValue           = uicontrol('Parent',hp,'String','0','style','edit',      'Position',[120 110 70  20],'backgroundcolor','white','Callback',@(src,evt)obj.set_probe_bias(src,evt));
			%% initialize parameters menu
			inp.UV_factor_text                   = uicontrol('Parent',hp,'String','UV factor',             'Position',[0   80 60 30]);
			inp.UV_factor_value                  = uicontrol('Parent',hp,'String',num2str(1),              'Position',[70  80 100 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
			inp.Rsun_text                        = uicontrol('Parent',hp,'String','Rsun [AU]',             'Position',[0   55 60 30]);
			inp.Rsun_value                       = uicontrol('Parent',hp,'String',num2str(1),              'Position',[70  55 100 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
			inp.U_text                           = uicontrol('Parent',hp,'String','U [V]',                 'Position',[0   30 60 30]);
			inp.U_value                          = uicontrol('Parent',hp,'String','',                      'Position',[70  30 100 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
			inp.update                           = uicontrol('Parent',hp,'String','Update',                'Position',[0   0 60 30],'Callback','lp.sweep_gui(''update'')');
			inp.reset                            = uicontrol('Parent',hp,'String','Reset',                 'Position',[70  0 60 30],'callback','lp.sweep_gui(''initialize'')');
			%% initialize s/c menu
			hsc = uipanel('Title','Spacecraft','FontSize',12,'BackgroundColor',[.95 1 1],'Position',[.7 .39 .3 .35]);
			inp.flag_sc                                    = uicontrol('Parent',hsc,'style','radio','String','Model spacecraft','Value',0,             'Position',[0   205 120 25]);
			inp.sc.example                                 = uicontrol('Parent',hsc,'String','Example spacecraft|Cluster|Solar Orbiter|THEMIS|Cassini','Position',[0   180 150 25],'style','popup','backgroundcolor','white','Callback', @setscexample);
			surf=lp.photocurrent;probtxt='spacecraft surface';for ii=1:numel(surf),probtxt(end+1:end+1+numel(surf{ii}))=['|' surf{ii}];end
			inp.sc.surface                                 = uicontrol('Parent',hsc,'String',probtxt,                                                  'Position',[0   155 150 30],'style','popup','backgroundcolor','white');
			inp.sc.total_area_text                         = uicontrol('Parent',hsc,'String','Total area [m2]',                                        'Position',[0   135 120 25]);
			inp.sc.total_area_value                        = uicontrol('Parent',hsc,'String',num2str(0),'style','edit',                 'Position',[120 135 50 25],'backgroundcolor','white');
			inp.sc.sunlit_area_text                        = uicontrol('Parent',hsc,'String','Sunlit area [m2]',                                       'Position',[0   110 120 25]);
			inp.sc.sunlit_area_value                       = uicontrol('Parent',hsc,'String',num2str(0),'style','edit',                'Position',[120 110 50 25],'backgroundcolor','white');
			inp.sc.antenna_guard_area_text                 = uicontrol('Parent',hsc,'String','Sunlit guard area [m2]',                                 'Position',[0   85 120 25],'Tooltipstring','Cross section area of pucks and guards, assuming similar photoelectron emission as antenna');
			inp.sc.antenna_guard_area_value                = uicontrol('Parent',hsc,'String','0','style','edit',                                       'Position',[120 85 50 25],'backgroundcolor','white');
			inp.sc.probe_refpot_as_fraction_of_scpot_text  = uicontrol('Parent',hsc,'String','Probe refpot/scpot',                                     'Position',[0   60 120 25],'Tooltipstring','The ratio between the probe reference potential and satellite potential');
			inp.sc.probe_refpot_as_fraction_of_scpot_value = uicontrol('Parent',hsc,'String',num2str(0),         'Position',[120 60 50 25],'style','edit','backgroundcolor','white');
			inp.sc.number_of_probes_text                   = uicontrol('Parent',hsc,'String','Number of probes',                                       'Position',[0   35 120 25]);
			inp.sc.number_of_probes_value                  = uicontrol('Parent',hsc,'String',num2str(0),                          'Position',[120 35 50 25],'style','edit','backgroundcolor','white');
			inp.sc.probe_distance_to_spacecraft_text       = uicontrol('Parent',hsc,'String','distance probe-sc [m]',                                  'Position',[0   10 120 25]);
			inp.sc.probe_distance_to_spacecraft_value      = uicontrol('Parent',hsc,'String',num2str(0),              'Position',[120 10 50 25],'style','edit','backgroundcolor','white');
			%% initialize plasma menu
			hpl= uipanel('Title','Plasma','FontSize',12,'BackgroundColor',[1 1 .95],'Position',[.7 .74 .3 .2]);
			inp.n         = uicontrol('Parent',hpl,'String','Ne [cc]',         'Position',[0 0 80 25]);
			inp.n_value   = uicontrol('Parent',hpl,'String','1','style','edit','Position',[80 0 90 25],'backgroundcolor','white');
			inp.T         = uicontrol('Parent',hpl,'String','T [eV]',          'Position',[0 25 80 25]);
			inp.T_value   = uicontrol('Parent',hpl,'String','1 1','style','edit','Position',[80 25 90 25],'backgroundcolor','white');
			inp.m         = uicontrol('Parent',hpl,'String','m [mp],0=me',     'Position',[0 50 80 25]);
			inp.m_value   = uicontrol('Parent',hpl,'String','0 1','style','edit','Position',[80 50 90 25],'backgroundcolor','white');
			inp.q         = uicontrol('Parent',hpl,'String','q [e]',           'Position',[0 75 80 25]);
			inp.q_value   = uicontrol('Parent',hpl,'String','-1 1','style','edit','Position',[80 75 90 25],'backgroundcolor','white');
			inp.vsc       = uicontrol('Parent',hpl,'String','Vsc [m/s]',       'Position',[0 100 80 25]);
			inp.vsc_value = uicontrol('Parent',hpl,'String','0','style','edit','Position',[80 100 90 25],'backgroundcolor','white');
			%% initialize plot menu
			hpl= uipanel('Title','Top panel','FontSize',12,'BackgroundColor',[1 1 .95],'Position',[.7 .94 .3 .06]);
			inp.toppanel.plot = uicontrol('Parent',hpl,'String','Resistance|Satellite IU|Antenna noise','Position',[0 0 150 25],'style','popup','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
			
			ud.inp=inp;
			obj.figHandle = figH;
			obj.UserData = ud;
			
		end
		function set_probe_type(obj,varargin)
			if nargin == 2, %set_probe_type(obj,probeId)
				probeId = varargin{1};
				set(obj.UserData.inp.probe.typeValue,'Value',probeId);
			elseif nargin == 3, %set_probe_type(obj,hEvent,event)
				hEvent = varargin{1};
				event = varargin{2};
				disp(event);
				idProbe = get(hEvent,'Value');
				obj.probeUsed = idProbe;
				probeParameters = obj.ProbeList(idProbe);
				indSurface = find(strcmp(probeParameters.surface,lp.photocurrent))+1;
				if indSurface,
					set(obj.UserData.inp.probe.surfaceValue,'Value',indSurface);
				else
					irf.log('critical',[surface '''' probeParameters.surface ''' is unknown by lp.photocurrent.']);
				end
				obj.set_probe_radius_sphere(obj.ProbeList(idProbe).radiusSphere);
				obj.set_probe_radius_wire(  obj.ProbeList(idProbe).radiusWire);
				obj.set_probe_length_wire(  obj.ProbeList(idProbe).lengthWire);
			end
		end
		function set_probe_radius_sphere(obj,radiusSphere)
			set(obj.UserData.inp.probe.radiusSphereValue,'String',num2str(radiusSphere*100,2)); % in cm
		end
		function get_probe_radius_sphere(obj)
			radiusSphereCm = get(obj.UserData.inp.probe.radiusSphereValue,'String'); % in cm
			radiusSphere = str2double(radiusSphereCm)*1e-2;
			if (obj.probeUsed ~= 1) && (radiusSphere ~= obj.ProbeList(obj.probeUsed).radiusSphere)
				obj.probeUsed = 1; % user defined
				obj.set_probe_type(1);
			end
			obj.ProbeList(obj.probeUsed).radiusSphere = radiusSphere;
		end
		function set_probe_radius_wire(obj,radiusWire)
			set(obj.UserData.inp.probe.radiusWireValue,'String',num2str(radiusWire*1e3,2)); % in mm
		end
		function set_probe_length_wire(obj,lengthWire) % input in [m]
			set(obj.UserData.inp.probe.lengthWireValue,'String',num2str(lengthWire*1e2,3)); % in [cm]
		end
		function set_probe_bias(obj,biasCurrent) % input in [A]
			set(obj.UserData.inp.probe.biasCurrentValue,'String',num2str(biasCurrent*1e-6,3)); % in [uA]
		end
	end
	methods (Static)
		function popupText = popup_list(inp)
			if iscell(inp) && (numel(inp) > 0)
				popupText =inp{1};
				for ii=2:numel(inp),
					popupText(end+1:end+1+numel(inp{ii}))=['|' inp{ii}];
				end
			elseif numel(inp) > 0 && all(isprop(inp,'name'))
				popupText ='';
				for ii=1:numel(inp),
					probeName = inp(ii).name;
					if ~isempty(probeName)
						popupText(end+1:end+1+numel(probeName))=[probeName '|'];
					end
				end
				popupText(end) =[];
			end
		end
	end
end

