function val = get(cdb,prop_name)
% GET Get ClusterDB properties from the specified object
% and return the value
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev

switch prop_name
case 'db'
	val = cdb.db;
case 'dp'
	val = cdb.dp;
case 'sp'
	val = cdb.sp;
otherwise
	error([prop_name,' Is not a valid ClusterDB property'])
end
