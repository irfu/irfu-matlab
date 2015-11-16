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

scList = {'Cluster',       @sc_Cluster;...
          'THEMIS',        @sc_THEMIS;...
          'Solar_Orbiter', @sc_Solar_Orbiter;...
          'THOR_SDP',      @sc_THOR_SDP;...
          'THOR_HFA',      @sc_THOR_HFA;...
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
		Sc.surface = 'solar cells';
		Sc.areaTotal = 25.66;
		Sc.areaSunlit = 3.87;
		Sc.areaSunlitGuard = 0.039;
		Sc.probeRefPotVsSatPot = 0.2;
		Sc.nProbes = 4;
		Sc.probeDistanceToSpacecraft = 44;
	end
	function Sc = sc_THEMIS
		Sc = lp.spacecraft;
		Sc.name  = 'THEMIS';
		Sc.probe = lp.default_lprobe('Cluster');
		Sc.surface = 'themis';
		Sc.areaTotal = NaN;
		Sc.areaSunlit = NaN;
		Sc.areaSunlitGuard = NaN;
		Sc.probeRefPotVsSatPot = 0.2;
		Sc.nProbes = 4;
		Sc.probeDistanceToSpacecraft = 30;
	end
	function Sc = sc_Solar_Orbiter
		Sc = lp.spacecraft;
		Sc.name  = 'Solar Orbiter';
		Sc.probe = lp.default_lprobe('Solar_Orbiter');
		Sc.surface = 'solar cells';
		Sc.areaTotal = 10;
		Sc.areaSunlit = 3;
		Sc.areaSunlitGuard = 0;
		Sc.probeRefPotVsSatPot = 0.2;
		Sc.nProbes = 3;
		Sc.probeDistanceToSpacecraft = 6;
	end
	function Sc = sc_THOR_SDP
		Sc = lp.spacecraft;
		Sc.name  = 'THOR SDP';
		Sc.probe = lp.default_lprobe('THOR_SDP');
		Sc.surface = 'themis';
		% 3.7m diameter, 0.7m height 
		Sc.areaSunlit = pi*(3.7/2)^2; % 10m^2
		Sc.areaTotal = 2*Sc.areaSunlit+pi*3.7*0.7; % 30m^2
		Sc.areaSunlitGuard = 0; % TO CHECK
		Sc.probeRefPotVsSatPot = 0.2;
		Sc.nProbes = 4;
		Sc.probeDistanceToSpacecraft = 50;
	end
	function Sc = sc_THOR_HFA
		Sc = lp.spacecraft;
		Sc.name  = 'THOR HFA';
		Sc.probe = lp.default_lprobe('THOR_HFA');
		Sc.surface = 'themis';
		Sc.areaSunlit = pi*(3.7/2)^2; % 10m^2
		Sc.areaTotal = 2*Sc.areaSunlit+pi*3.7*0.7; % 30m^2
		Sc.areaSunlitGuard = 0; % TO CHECK
		Sc.probeRefPotVsSatPot = 0.2;
		Sc.nProbes = 6;
		Sc.probeDistanceToSpacecraft = 1;
	end
end

