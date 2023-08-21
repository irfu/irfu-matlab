function dir_list = caa_get_subdirs(iso_t, dt, cl_id)
%CAA_GET_SUBDIRS     Get the data subdirectories for all subintervals during
%                    the given interval, and for the given Cluster spacecraft.
%
% res = caa_get_subdirs(st, dt, cl_id)
% res = caa_get_subdirs(iso_t, dt, cl_id)

% Copyright 2008 Mikael Lundberg

narginchk(3,3);

if ~isnumeric(cl_id) || cl_id < 1 || cl_id > 4
  error('Wrong Cluster id. Valid values are: [1-4]')
end

SPLIT_INT = 3; % 3 hour subintervals
if ismac
  BASE_DIR = '/Volumes/caa/l1';
else
  BASE_DIR = '/data/caa/l1';
end

[st,dt] = irf_stdt(iso_t,dt);

t = fromepoch(st);
t0 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);
t = fromepoch(st+dt);
t1 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);
if t1>=st+dt, t1 = t1 - SPLIT_INT*3600; end

old_pwd = pwd;
dir_list = {};
good_dir = {};
for t=t0:SPLIT_INT*3600:t1
  y = fromepoch(t);
  main_int = [BASE_DIR '/' num2str(y(1)) '/' irf_fname(t) '/C' num2str(cl_id)];
  if ~exist(main_int,'dir'), continue, end

  cd(main_int)
  d = dir('*_*');
  if isempty(d), continue, end
  %	good_dir = {};
  for j=1:length(d)
    if ~d(j).isdir, continue, end
    if caa_is_valid_dirname(d(j).name), good_dir{end+1} = [main_int '/' d(j).name]; end
  end
  %	if isempty(good_dir), continue, end
end
cd(old_pwd)
dir_list = sort(good_dir);

%	for j=1:length(good_dir)
%		subdir = [main_int '/' good_dir{j}];
%		cd(subdir)
%		[st_t,dt_tmp] = caa_read_interval();
%		if isempty(st_t), continue, end









%dir_list = [];
%BASE_DIR = '/data/caa/l1';
%
%if caa_is_valid_dirname(main_int)
%   main_int = [BASE_DIR '/' main_int(1:4) '/' main_int];
%end
%
%if ~exist(main_int, 'dir'), return, end
%
%good_dir = {};
%old_pwd=pwd;
%
%cd(main_int)
%c = dir('C*');
%if isempty(c), error(['Found no valid subdirs for: ' main_int]), end
%for j = 1:length(c)
%   if ~c(j).isdir, continue, end
%
%   d = dir([c(j).name '/*_*']);
%	if isempty(d), continue, end
%
%	for k=1:length(d)
%		if ~d(k).isdir, continue, end
%		if caa_is_valid_dirname(d(k).name)
%		   good_dir = {good_dir{:} [c(j).name '/' d(k).name]};
%		end
%	end
%end
%cd(old_pwd)
%dir_list = good_dir;