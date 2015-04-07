function plasmaList = default_plasma( plasmaNames )
%LP.DEFAULT_PLASMA different default plasma environments

plasmaNamesList = {'sw@1AU','plasma sheet'};
PlasmaFunc      = {@p_sw_1AU, @p_plasma_sheet,};

if nargin == 0 && nargout == 0
	for iPlasma = 1:numel(plasmaNamesList),
		disp([num2str(iPlasma) '. ' plasmaNamesList{iPlasma}]);
	end
	return;
elseif nargin == 0
	plasmaNames = plasmaNamesList; 
end
if ischar(plasmaNames), plasmaNames = {plasmaNames};   end

iFoundPlasma = [];
for iPlasma = 1:numel(plasmaNames),
	iFoundPlasma = [iFoundPlasma find(strcmp(plasmaNames(iPlasma),plasmaNamesList))]; %#ok<AGROW>
end

if iFoundPlasma
	plasmaList = PlasmaFunc{iFoundPlasma(1)}();
	for j=2:numel(iFoundPlasma)
		plasmaList(j) = PlasmaFunc{iFoundPlasma(j)}();
	end
end

% 	properties
% 		q
% 		n
% 		mp
% 		T
% 		v
% 	end

	function Plasma = p_sw_1AU
		Plasma = lp.plasma;
		Plasma.name = 'SW @1AU';
		Plasma.qe = [-1 1];
		Plasma.n = 4e6;
		Plasma.mp = [0 1];
		Plasma.T = 10;
		Plasma.v = 400e3;
	end
	function Plasma = p_plasma_sheet
		Plasma = lp.plasma;
		Plasma.name = 'Plasma sheet';
		Plasma.qe = [-1 1];
		Plasma.n = 1e6;
		Plasma.mp = [0 1];
		Plasma.T = [1e3 5e3];
		Plasma.v = 0;
	end
end

