function ScListOutput = default_spacecraft(scNames)
%LP.DEFAULT_SPACECRAFT return existing default spacecraft
%
% Spacecraft parameters are defined by LP.SPACECRAFT
%
%  LP.SPACECRAFT print the list of available spacecraft
%
%  OUT = LP.SPACECRAFT return all cell array with all spacecraft
%
%  OUT = LP.SPACECRAFT({name1,name2,..}) return selected spacecraft
%

scList = {'Cluster',       @sc_Cluster,...
          'THEMIS',        @sc_THEMIS;...
          'Solar_Orbiter', @sc_Solar_Orbiter;...
					};
scNamesList = scList(:,1);
ScFunc      = scList(:,2);

if nargin == 0 && nargout == 0
	for iSc = 1:numel(scNamesList),
		disp([num2str(iSc) '. ' scNamesList{iSc}]);
	end
	return;
elseif nargin == 0
	scNames = scNamesList;
end
if ischar(scNames), scNames = {scNames};   end

iFoundSc = [];
for iSc = 1:numel(scNames),
	iFoundSc = [iFoundSc find(strcmp(scNames(iSc),scNamesList))]; %#ok<AGROW>
end

if iFoundSc
	ScListOutput = ScFunc{iFoundSc(1)}();
	for j=2:numel(iFoundSc)
		ScListOutput(j) = ScFunc{iFoundSc(j)}();
	end
end

	function Sc = sc_Cluster
		Sc = lp.spacecraft;
		Sc.name  = 'Cluster';
		Sc.probe = lp.default_lprobe('Cluster');
		Sc.surface = 'solar cells';%    set(data.inp.sc.surface,'Value',find(strcmp('solar cells',lp.photocurrent))+1); % solar cells
		Sc.areaTotal = 25.66; %   set(data.inp.sc.total_area_value,'string','25.66');
		Sc.areaSunlit = 3.87; %   set(data.inp.sc.sunlit_area_value,'string','3.87');
		Sc.areaSunlitGuard = 0.039; %    set(data.inp.sc.antenna_guard_area_value,'string','0.039');
		Sc.probeRefPotvsSatPot = 0.2;		%     set(data.inp.sc.probe_refpot_as_fraction_of_scpot_value,'string','.2');
		Sc.nProbes = 4; %     set(data.inp.sc.number_of_probes_value,'string','4');
		Sc.probeDistanceToSpacecraft = 44; %     set(data.inp.sc.probe_distance_to_spacecraft_value,'string','44');
	end
	function Sc = sc_THEMIS
		Sc = lp.spacecraft;
		Sc.probe = lp.default_lprobe('Cluster'); %data.probe.type='spherical';
	end
	function Sc = sc_Solar_Orbiter
		Sc = lp.spacecraft;
		Sc.name  = 'Solar Orbiter';
		Sc.probe = lp.default_lprobe('Solar_Orbiter');
		Sc.surface = 'solar cells';% 
		Sc.areaTotal = 10; %  
		Sc.areaSunlit = 3; % 
		Sc.areaSunlitGuard = 0; 
		Sc.probeRefPotvsSatPot = 0.2;	
		Sc.nProbes = 3; %    
		Sc.probeDistanceToSpacecraft = 6;
	end
end

