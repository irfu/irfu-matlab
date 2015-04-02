function ui( varargin )
%LP.UI User interface to +lp package
%   LP.UI(spacecraft,probe,plasma,parameters)

Units=irf_units;
persistent message;
if isempty(message), % run only the first time during the session
    message='You can anytime access all the results from get(gcf,''userdata'').';
    disp(message);
end

%% Defaults
SpacecraftList = lp.default_spacecraft;
ProbeList      = lp.default_probes;
PlasmaList     = lp.default_plasma_models;

%% Check input
while ~isempty(args)
	if isa(args{1},'spacecraft'),
		SpacecraftList = args{1};
		args(1) =[];
	elseif isa(args{1},'lpProbe'),
		ProbeList = args{1};
		args(1) =[];
	elseif isa(args{1},'plasma')
		PlasmaList = args{1};
		args(1) =[];
	else
		
	end
end
if      nargin == 0, action='initialize';end

end

