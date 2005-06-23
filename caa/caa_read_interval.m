function [st,dt] = caa_read_interval(sp)
%CAA_READ_INTERVAL  read caa file .interval
%
% [st,dt] = caa_read_interval([sp])
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

if nargin<1, sp=pwd; end

old_pwd = pwd;
cd(sp);
if exist('./.interval','file')
	[st_s,dt_s] = textread('./.interval','%s %s',-1);
	st_s = st_s{1}; dt_s = dt_s{1};
	if nargout==0
		disp([st_s ' : ' dt_s ' sec'])
	elseif nargout==1
		st = st_s;
	elseif nargout==2
		dt = str2num(dt_s);
		st = st_s;
	else
		error('wrong number of input arguments')
	end
else
	if nargout==0
		disp('cannot find .version')
	elseif nargout==1
		st = '';
	elseif nargout==2
		dt = [];
		st = '';
	else
		error('wrong number of input arguments')
	end
end
cd(old_pwd)
