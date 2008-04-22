function display(dobj,mode)
%DISPLAY(dobj)  display a DATAOBJ object
%
% DISPLAY(dobj,[mode])
%
% mode : 'full'  - list everything
%        'short' - only data varibles
%
% $Revision$  $Date$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin < 2, m = 0;
else
	switch lower(mode)
		case {'f','full'}
			m = 1;
		case {'s','short'}
			m = 0;
		otherwise
			error('invalid value for MODE')
	end
end

disp(' ')
disp(['dataobj object created : ' dobj.FileModDate])
disp(' ')
disp('Variables:')
nvars = size(dobj.vars,1);
if nvars>0
	for v=1:nvars
		if m == 0 && strcmpi(dobj.data.(dobj.vars{v,1}).type,'char'), continue, end
		disp([dobj.vars{v,1} ' : ' dobj.data.(dobj.vars{v,1}).type ' : '...
			num2str(dobj.data.(dobj.vars{v,1}).nrec) ' recs' ])
	end
else
	disp('empty')
end
disp(' ')	
