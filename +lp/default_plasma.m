function plasmaListOutput = default_plasma( plasmaNames )
%LP.DEFAULT_PLASMA different default plasma environments
%
% Plasma models are defined by class LP.PLASMA
%
%  LP.DEFAULT_PLASMA print the list of available plasma models
%
%  OUT = LP.DEFAULT_PLASMA return cell array with all plasma models
%
%  OUT = LP.DEFAULT_PLASMA({name1,name2,..}) return selected plasma models
%

plasmaList = {'sw@1AU',       @p_sw_1AU;...
              'swfast@1AU',   @p_sw_fast_1AU;...
              'swslow@1AU',   @p_sw_slow_1AU;...
	            'sw@028AU',     @p_sw_028AU;...
	            'plasma sheet', @p_plasma_sheet...
             };
plasmaNamesList = plasmaList(:,1);
PlasmaFunc      = plasmaList(:,2);

if nargin == 0 && nargout == 0
	for iPlasma = 1:numel(plasmaNamesList)
		disp([num2str(iPlasma) '. ' plasmaNamesList{iPlasma}]);
	end
	return;
elseif nargin == 0
	plasmaNames = plasmaNamesList; 
end
if ischar(plasmaNames), plasmaNames = {plasmaNames};   end

iFoundPlasma = [];
for iPlasma = 1:numel(plasmaNames)
	iFoundPlasma = [iFoundPlasma find(strcmp(plasmaNames(iPlasma),plasmaNamesList))]; %#ok<AGROW>
end

if iFoundPlasma
	plasmaListOutput = PlasmaFunc{iFoundPlasma(1)}();
	for j=2:numel(iFoundPlasma)
		plasmaListOutput(j) = PlasmaFunc{iFoundPlasma(j)}();
	end
end

% 	properties
% 		q  [e]
% 		n  [m^-3]
% 		mp [mp] if 0 then e-
% 		T  [eV]
% 		v  [m/s]
% 	end

	function Plasma = p_sw_1AU
		Plasma      = lp.plasma;
		Plasma.name = 'SW @1AU';
		Plasma.qe   = [-1 1];
		Plasma.n    = 4e6;
		Plasma.mp   = [0 1];
		Plasma.T    = 10;
		Plasma.v    = 400e3;
	end
	function Plasma = p_sw_fast_1AU
		Plasma      = lp.plasma;
		Plasma.name = 'SW fast @1AU';
		Plasma.qe   = [-1 1];
		Plasma.n    = 3e6;
		Plasma.mp   = [0 1];
		Plasma.T    = 20;
		Plasma.v    = 600e3;
	end
	function Plasma = p_sw_slow_1AU
		Plasma      = lp.plasma;
		Plasma.name = 'SW slow @1AU';
		Plasma.qe   = [-1 1];
		Plasma.n    = 10e6;
		Plasma.mp   = [0 1];
		Plasma.T    = 4;
		Plasma.v    = 400e3;
	end
	function Plasma = p_sw_028AU
		Plasma      = lp.plasma;
		Plasma.name = 'SW+strahl @0.28AU';
		Plasma.qe   = [-1  -1  1];
		Plasma.n    = [100 8  108]*1e6;
		Plasma.mp   = [0   0   1];
		Plasma.T    = [25 100 30];
		Plasma.v    = 400e3;
	end
	function Plasma = p_plasma_sheet
		Plasma      = lp.plasma;
		Plasma.name = 'Plasma sheet';
		Plasma.qe   = [-1 1];
		Plasma.n    = 1e6;
		Plasma.mp   = [0 1];
		Plasma.T    = [1e3 5e3];
		Plasma.v    = 0;
	end
end

