function val = get(O,prop_name)
% GET Get ClusterProc properties from the specified object
% and return the value
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev

switch prop_name
case 'sp'
	val = O.sp;
otherwise
	error([prop_name,' Is not a valid ClusterProc property'])
end
