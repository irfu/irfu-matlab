function ScList = default_spacecraft(scNames)
%LP.DEFAULT_SPACECRAFT return existing default spacecraft
scNamesList = {'Cluster','THEMIS'};
ScFunc      = {@sc_Cluster @sc_THEMIS};

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
	iFoundSc = [iFoundSc find(strcmp(scNames(iSc),scNamesList))];
end

if iFoundSc
	ScList = ScFunc{iFoundSc(1)}();
	for j=2:numel(iFoundSc)
		ScList(j) = ScFunc{iFoundSc(j)}();
	end
end

	function Sc = sc_Cluster
		Sc = lp.spacecraft;
		Sc.name  = 'Cluster';
		Sc.probe = 'spherical'; %data.probe.type='spherical';
		%    set(data.inp.probe.type,'Value',1);
		%    data.probe.surface='themis';
		%    set(data.inp.probe.surface,'Value',find(strcmp('cluster',lp.photocurrent))+1);
		%    set(data.inp.probe.length_value,'style','text','string','');
		%    set(data.inp.probe.radius_value,'string','4');
		Sc.surface = 'solar cells';%    set(data.inp.sc.surface,'Value',find(strcmp('solar cells',lp.photocurrent))+1); % solar cells
		Sc.areaTotal = 25.66; %   set(data.inp.sc.total_area_value,'string','25.66');
		Sc.areaSunlit = 3.87; %   set(data.inp.sc.sunlit_area_value,'string','3.87');
		Sc.areaSunlitGuard = 0.039; %    set(data.inp.sc.antenna_guard_area_value,'string','0.039');
		%     set(data.inp.sc.probe_refpot_as_fraction_of_scpot_value,'string','.2');
		%     set(data.inp.sc.number_of_probes_value,'string','4');
		%     set(data.inp.sc.probe_distance_to_spacecraft_value,'string','44');
		%     set(data.inp.Rsun_value,'string','1');
		%     set(data.inp.probe.total_vs_sunlit_area_value,'string','4');
		%     data.probe.total_vs_sunlit_area=4;
		%     set(data.inp.n_value,'string','1');
		%     set(data.inp.T_value,'string','100 500');
	end
	function Sc = sc_THEMIS
		Sc = lp.spacecraft;
		Sc.probe = 'spherical'; %data.probe.type='spherical';
	end
end

