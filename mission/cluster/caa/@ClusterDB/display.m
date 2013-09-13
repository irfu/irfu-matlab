function display(cdb)
% DISPLAY(cdb) Display a ClusterDB object
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev

stg = sprintf(...
'ClusterDB:\ndb-> %s\ndp-> %s\nsp-> %s', ...
cdb.db, cdb.dp, cdb.sp);
disp(stg)
