function [data, ok, msg] = caa_get(iso_t,dt,cl_id,var_name,ops_s,varargin)
%CAA_GET  read data from caa Matlab files
%
% data = caa_get(start_t,dt/stop_t,cl_id,var_name,ops_s)
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if ismac
  DP = '/Volumes/caa/l1';    % TODO: Add option to change base dir via input param! (ML)
else
  DP = '/data/caa/l1';
end
data = [];
ok=0;
msg='';

NEED_SAME_TM = 0;
DO_MEAN = 0;
HAVE_LOAD_ARGS = 0;

if nargin > 4
  if nargin > 5
    args = varargin;
  end
  switch ops_s
    case 'need_same_tm'
      NEED_SAME_TM = 1;
    case 'mean'
      DO_MEAN = 1;
    case 'base_dir'
      DP = args{1};
    case 'load_args'
      HAVE_LOAD_ARGS = 1;
      load_args = args;
      eval_str = 'var_name,cl_id';
      for k=1:numel(load_args)
        eval_str = [eval_str ',load_args{' num2str(k) '}']; %#ok<AGROW>
      end
    otherwise
      irf_log('fcal','unknown option')
  end
end

SPLIT_INT = 3; % 3 hour subintervals

[st,dt] = irf_stdt(iso_t,dt);

t = fromepoch(st);
t0 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);
t = fromepoch(st+dt);
t1 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);
if t1>=st+dt, t1 = t1 - SPLIT_INT*3600; end

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
    if caa_is_valid_dirname(d(j).name), good_dir = [good_dir {d(j).name}]; end %#ok<AGROW>
  end
  if isempty(good_dir), continue, end
  
  for j=1:length(good_dir)
    subdir = [main_int '/' good_dir{j}];
    cd(subdir)
    [st_t,dt_tmp] = caa_read_interval();
    if isempty(st_t), continue, end
    st_tmp = iso2epoch(st_t);
    if (st_tmp+dt_tmp <= st) || (st_tmp >= st+dt), continue, end % subinterval starts before/after the interval
    [ok,tm] = c_load('mTMode?',cl_id);
    if ~ok, continue, end
    ttt.st = st_tmp;
    ttt.dt = dt_tmp;
    
    if NEED_SAME_TM && tm~=tm(1)*ones(size(tm))
      error('tape mode changes during the selected time inteval')
    end
    ttt.mode = tm(1);
    ttt.dir = subdir;
    if isempty(mode_list), mode_list = ttt;
    else, mode_list(end+1) = ttt; %#ok<AGROW>
    end
  end
  clear good_dir
end
if isempty(mode_list), cd(old_pwd), return, end

% Concatenate intervals
[starts,ii] = sort([mode_list.st]);
dts = [mode_list.dt]; dts = dts(ii);
for j = ii
  if regexp(var_name,'^mTMode([1-4]|?)$')==1
    tt = [mode_list(j).st mode_list(j).mode];
  else
    cd(mode_list(j).dir);
    if HAVE_LOAD_ARGS
      eval(['[ok, tt, msg] = c_load(' eval_str ');']);
    else
      [ok, tt, msg] = c_load(var_name,cl_id);
    end
    if ~ok || isempty(tt), continue, end
    % Remove NaN times
    % TODO: times must never be NaN.
    if ~isstruct(tt), tt(isnan(tt(:,1)),:) = []; end
    if isempty(tt), continue, end
    
    % Uncorrect delta offsets. This is black magic...
    if regexp(var_name,'^diEs([1-4]|?)p(12|32|34)$')==1
      % Delta offsets
      pp = str2double(var_name(end-1:end));
      [ok,Delauto] = c_load('D?p12p34',cl_id);
      if ~ok || isempty(Delauto)
        irf_log('load',irf_ssub('Cannot load/empty D?p12p34',cl_id))
      else
        tt = caa_corof_delta(tt,pp,Delauto,'undo');
      end
    end
    
    % Append time to variables which does not have it
    % 946684800 = toepoch([2000 01 01 00 00 00])
    if ~isstruct(tt) && ( tt(1,1) < 946684800 )
      if size(tt,1)==1, tt = [starts(j)+dts(j)/2.0 tt]; %#ok<AGROW>
      else, error('loaded bogus data')
      end
    elseif DO_MEAN
      mm = mean(tt(~isnan(tt(:,2)),2:end));
      if any(~isnan(mm)), tt = [starts(j)+dts(j)/2.0 mm];
      else, continue
      end
    end
  end
  
  if isempty(data), data = tt;
  elseif ~(isstruct(data) || isstruct(tt))
    data = caa_append_data(data,tt);
  elseif (isstruct(data) && isstruct(tt))
    if ~isfield(data, 'int1')
      temp.int1 = data;
      data = temp;
      data.int1.mode_list = mode_list(1);
      clear temp;
    end
    data.(['int' num2str(j)]) = tt;
    data.(['int' num2str(j)]).mode_list = mode_list(j);
  end
end

if nargout > 1, ok = ~isempty(data); end

cd(old_pwd)
