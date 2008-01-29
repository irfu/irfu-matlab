function data = caa_get(iso_t,dt,cl_id,var_name,ops_s)
%CAA_GET  read data from caa Matlab files
%
% data = caa_get(start_t,dt/stop_t,cl_id,var_name,ops_s)
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

DP = '/data/caa/l1';
data = [];

NEED_SAME_TM = 0;

if nargin > 4
	if strcmp(ops_s, 'need_same_tm'), NEED_SAME_TM = 1; end
end

SPLIT_INT = 3; % 3 hour subintervals

[st,dt] = irf_stdt(iso_t,dt);

t = fromepoch(st);
t0 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);
t = fromepoch(st+dt);
t1 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);

old_pwd = pwd;
mode_list = [];
for t=t0:SPLIT_INT*3600:t1
	y = fromepoch(t);
	main_int = [DP '/' num2str(y(1)) '/' irf_fname(t) '/C' num2str(cl_id)];
	if ~exist(main_int,'dir'), continue, end
	
	%disp(main_int)
	cd(main_int)
	d = dir('*_*');
	if isempty(d), continue, end
	good_dir = {};
	for j=1:length(d)
		if ~d(j).isdir, continue, end
		if caa_is_valid_dirname(d(j).name), good_dir = {good_dir{:} d(j).name}; end
	end
	if isempty(good_dir), continue, end
	
	for j=1:length(good_dir)
		subdir = [main_int '/' good_dir{j}];
		cd(subdir)
		[st_t,dt_tmp] = caa_read_interval();
		[ok,tm] = c_load('mTMode?',cl_id);
		if isempty(st_t) || ~ok, continue, end
		ttt.st = iso2epoch(st_t);
		ttt.dt = dt_tmp;

		if NEED_SAME_TM && tm~=tm(1)*ones(size(tm))
			error('tape mode changes during the selected time inteval')
		end
		ttt.mode = tm(1);
		ttt.dir = subdir;
		if isempty(mode_list), mode_list = ttt;
        else mode_list(end+1) = ttt; %#ok<AGROW>
		end
	end
	clear good_dir
end
if isempty(mode_list), cd(old_pwd), return, end

% Concatenate intervals
[starts,ii] = sort([mode_list.st]);
for j = ii;
	cd(mode_list(j).dir);
	[ok, tt] = c_load(var_name,cl_id);
	if ~ok || isempty(tt), continue, end
	% Remove NaN times
	% TODO: times must never be NaN.
	tt(isnan(tt(:,1)),:) = []; if isempty(tt), continue, end
	
	% Append time to variables which does not have it
	% 946684800 = toepoch([2000 01 01 00 00 00])
	if tt(1,1) < 946684800, tt = [starts(j) tt]; end %#ok<AGROW>
	
	if isempty(data), data = tt;
    else data = caa_append_data(data,tt);
	end
end

if ~isempty(data), data = irf_tlim(data,st +[0 dt]); end

cd(old_pwd)
