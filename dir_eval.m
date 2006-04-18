function dir_eval(eval_s,dirs)
%DIR_EVAL  evaluate in subdirectories
%
% DIR_EVAL(EVAL_S,[DIRS])
%    Run eval(EVAL_S) in all subdirectories of the curret directory or 
%    in the ones specified by DIRS (cell array)
%
% See also EVAL, C_EVAL
%
% $Id$

% Copyright 2006 Yuri Khotyaintsev

error(nargchk(1,2,nargin))

if nargin<2
	dirs = {};
	d = dir;
	for j=1:length(d)
		if ~strcmp(d(j).name,'.') && ~strcmp(d(j).name,'..')
			if isempty(dirs) dirs = { d(j).name };
			else dirs = [ dirs {d(j).name} ];
			end
		end
	end
end

if isempty(dirs), irf_log('fcal','Empty directory list'), return, end

old_pwd = pwd;
for j=1:length(dirs)
	lasterror('reset')
	try
		disp(['Executing in ' dirs{j}])
		cd (dirs{j})
		eval(eval_s)
	catch
		disp(['Error in ' dirs{j}])
		e = lasterror;
		disp(e.message)
	end
	cd (old_pwd)
end
