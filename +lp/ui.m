function ui( varargin )
%LP.UI User interface to +lp package
%   LP.UI(spacecraft,probe,plasma,parameters)

%% Message shown once during the session
Units=irf_units;
persistent message;
if isempty(message), % run only the first time during the session
    message='You can anytime access all the results from get(gcf,''userdata'').';
    disp(message);
end

%% Defaults
SpacecraftList = lp.default_spacecraft;
ProbeList      = lp.default_lprobe;
PlasmaList     = lp.default_plasma;

%% Check input
args=varargin;
while ~isempty(args)
	if isa(args{1},'lp.spacecraft'),
		SpacecraftList = args{1};
		args(1) =[];
	elseif isa(args{1},'lp.lprobe'),
		ProbeList = args{1};
		args(1) =[];
	elseif isa(args{1},'lp.plasma')
		PlasmaList = args{1};
		args(1) =[];
	else
		
	end
end

%% Initialize IDE
lp.gui_initialize(SpacecraftList, ProbeList, PlasmaList);

end

