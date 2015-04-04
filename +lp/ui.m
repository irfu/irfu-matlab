classdef ui
	%LP.UI User interface to +lp package
	%   LP.UI(spacecraft,probe,plasma,parameters)
	properties
		SpacecraftList 
		ProbeList      
		PlasmaList
		figHandle
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
		
		
	end
end
	
