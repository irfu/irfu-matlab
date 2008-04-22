function res = getv(dobj,var_s)
%GETV(dobj, var_s)  get a variable
%
% $Revision$  $Date$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,2,nargin))

if ~ischar(var_s), error('VAR_S must be a stirng'), end

% Take care of long variables (>63 symbols)
if length(var_s)>63
	var_s = [var_s(1:60) '...' var_s(end)];
	disp('trying truncated variable name')
end

nvars = size(dobj.vars,1);
if nvars>0
	for v=1:nvars
		if strcmpi(dobj.vars{v,1},var_s) || strcmpi(dobj.vars{v,2},var_s)
			res = dobj.data.(dobj.vars{v,1});
			res.name = var_s;
			return
		end
	end
end

error(['No such variable : ' var_s])

