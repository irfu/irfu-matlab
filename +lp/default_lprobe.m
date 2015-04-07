function lprobeList = default_lprobe( lprobeNames )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

lprobeNamesList = {'sphere',   'wire',  'sphere+wire',  'Cluster', 'MMS_SDP', 'MMS_ADP'};
LpFunc          = {@lp_sphere @lp_wire @lp_sphere_wire @lp_Cluster @lp_MMS_SDP @lp_MMS_ADP};

if nargin == 0 && nargout == 0
	for iLprobe = 1:numel(lprobeNamesList),
		disp([num2str(iLprobe) '. ' lprobeNamesList{iLprobe}]);
	end
	return;
elseif nargin == 0
	lprobeNames = lprobeNamesList; 
end
if ischar(lprobeNames), lprobeNames = {lprobeNames};   end

iFoundLprobe = [];
for iLprobe = 1:numel(lprobeNames),
	iFoundLprobe = [iFoundLprobe find(strcmp(lprobeNames(iLprobe),lprobeNamesList))]; %#ok<AGROW>
end

if iFoundLprobe
	lprobeList = LpFunc{iFoundLprobe(1)}();
	for j=2:numel(iFoundLprobe)
		lprobeList(j) = LpFunc{iFoundLprobe(j)}();
	end
end

% 		name
% 		surface
% 		radiusSphere
% 		radiusWire
% 		lengthWire

	function Lprobe = lp_sphere
		Lprobe = lp.lprobe;
		Lprobe.name = 'sphere';
		Lprobe.surface = 'themis';
		Lprobe.radiusSphere = 0.04; % 4cm
	end
	function Lprobe = lp_wire
		Lprobe = lp.lprobe;
		Lprobe.name = 'cylinder/wire';
		Lprobe.surface = 'themis';
		Lprobe.radiusWire = 1e-3;
		Lprobe.lengthWire = 1;
	end
	function Lprobe = lp_sphere_wire
		Lprobe = lp.lprobe;
		Lprobe.name = 'sphere+wire';
		Lprobe.surface = 'themis';
		Lprobe.radiusSphere = 0.04; % 4cm
		Lprobe.radiusWire = 1e-3;
		Lprobe.lengthWire = 1;
	end
	function Lprobe = lp_Cluster
		Lprobe = lp.lprobe;
		Lprobe.name = 'Cluster';
		Lprobe.surface = 'cluster';
		Lprobe.radiusSphere = 0.04; % 4cm
		Lprobe.radiusWire = 1e-3;
		Lprobe.lengthWire = 1;
	end
	function Lprobe = lp_MMS_SDP
		Lprobe = lp.lprobe;
		Lprobe.name = 'MMS SDP';
		Lprobe.surface = 'TiN';
		Lprobe.radiusSphere = 0.04; % 4cm
		Lprobe.radiusWire = 0.12e-3;
		Lprobe.lengthWire = 1.75;
	end
	function Lprobe = lp_MMS_ADP
		Lprobe = lp.lprobe;
		Lprobe.name = 'MMS ADP';
		Lprobe.surface = 'cluster';
		Lprobe.radiusWire = 0.5e-2;
		Lprobe.lengthWire = 1;
	end

end

