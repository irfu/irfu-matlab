classdef ui < handle
	%LP.UI User interface to +lp package
	%   LP.UI(spacecraft,probe,plasma,parameters)
	properties (SetAccess = protected)
		SpacecraftList
		SpacecraftUserDefined
		ProbeList
		ProbeUserDefined
		PlasmaList
		PlasmaUserDefined
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
			%% Initialize IDE
			obj.new_ide();
			obj.SpacecraftUserDefined = obj.SpacecraftList(1);
			obj.ProbeUserDefined = obj.ProbeList(1);
			obj.PlasmaUserDefined = obj.PlasmaList(1);
		end
		function set.PlasmaUserDefined(obj,PlasmaModel)
			obj.PlasmaUserDefined = PlasmaModel;
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
			hp = uipanel('Title','Probe','FontSize',12,'BackgroundColor',[1 0.95 1],'Position',[.7 .0 .3 .39]);
			popuptxt = obj.popup_list(obj.ProbeList);
			inp.probe.type                       = uicontrol('Parent',hp,'String',popuptxt,'style','popup','Position',[2 230 150 30],'backgroundcolor','white','Callback',@(src,ev)obj.setprobetype(src,ev));
			inp.probe.surface                    = uicontrol('Parent',hp,'String',lp.photocurrent,         'Position',[2 210 150 30],'style','popup','backgroundcolor','white');
			inp.probe.total_vs_sunlit_area_text  = uicontrol('Parent',hp,'String','total/sunlit area',     'Position',[0   185 120 30]);
			inp.probe.total_vs_sunlit_area_value = uicontrol('Parent',hp,'String','','style','text',       'Position',[120 185 70 30],'backgroundcolor','white');
			inp.probe.length_text                = uicontrol('Parent',hp,'String','probe length [cm]',     'Position',[0   160 120 30]);
			inp.probe.length_value               = uicontrol('Parent',hp,'String','','style','text',       'Position',[120 160 70 30],'backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
			inp.probe.radius_text                = uicontrol('Parent',hp,'String','probe radius [cm]',     'Position',[0   135 120 30]);
			inp.probe.radius_value               = uicontrol('Parent',hp,'String','',                      'Position',[120 135 70 30],'style','edit','backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
			inp.probe.bias_current_text          = uicontrol('Parent',hp,'String','bias current [uA]',     'Position',[0   110 120 30]);
			inp.probe.bias_current_value         = uicontrol('Parent',hp,'String','0','style','edit',      'Position',[120 110 70 30],'backgroundcolor','white','Callback','lp.sweep_gui(''update'')');
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
		function setprobetype(obj,hEvent,event) 
			disp(event);
			idProbe = get(hEvent,'Value');
			probeParameters = obj.ProbeList(idProbe);
			indSurface = find(strcmp(probeParameters.surface,lp.photocurrent));
			if indSurface,
				set(obj.UserData.inp.probe.surface,'Value',indSurface);
			else
				irf.log('critical',[surface '''' probeParameters.surface ''' is unknown by lp.photocurrent.']);
			end
%		radiusSphere
%		radiusWire     % if radius wire is two numbers than use stazer formula
%		lengthWire

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

