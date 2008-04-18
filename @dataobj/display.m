function display(dobj)
% DISPLAY(dobj) Display a dataobj object
%
% $Revision$  $Date$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


disp(' ')
disp(['dataobj object created : ' dobj.FileModDate])
disp(' ')
disp('Variables:')
nvars = size(dobj.vars,1);
if nvars>0
	for v=1:nvars
		disp([dobj.vars{v,1} ' : ' dobj.data.(dobj.vars{v,1}).type ' : '...
			num2str(dobj.data.(dobj.vars{v,1}).nrec)])
	end
else
	disp('empty')
end
disp(' ')	
