function O = set(O,varargin)
% SET Set ClusterProc properties and return the updated object
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev

property_argin = varargin;
while length(property_argin) >= 2
	prop = property_argin{1};
	val = property_argin{2};
	property_argin = property_argin(3:end);
	switch prop
	case 'sp'
		O.sp = val;
	otherwise
		error('ClusterProc properties: db, dp, sp')
	end
end
