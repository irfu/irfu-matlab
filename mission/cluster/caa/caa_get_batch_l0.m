function caa_get_batch_l0(iso_t,dt,cl_id,sdir,srcvars)
%CAA_GET_BATCH_L0  CAA get L0 data by CaaRunner
%
% caa_get_batch_l0(iso_t,dt,cl_id,sdir,srcvars)
%
% Input: iso_t - ISO epoch
%           dt - length of time interval in sec
%        cl_id - Cluster ID
%         sdir - directory to save the data
%      srcvars - see help ClusterProc/getData, ex. 'dies|die'
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(5,5)
sc_list = cl_id;

REQ_INT = 60; % Intervals (sec) for which we request FDM
SPLIT_INT = 90*60; % Typical interval length (sec)
MAX_SKIP = 2; % Number of NM frames we skip after BM interval

warning off  'ISDAT:serverWarning'
warning off  'ISDAT:serverMessage'

DB_S = c_ctl(0,'isdat_db');
DP_S = c_ctl(0,'data_path');

% Create the storage directory if it does not exist
if ~exist(sdir, 'dir')
  [SUCCESS,MESSAGE] = mkdir(sdir);
  if SUCCESS, irf_log('save',['Created storage directory ' sdir])
  else, error(MESSAGE)
  end
end

st = iso2epoch(iso_t);

% First we check if we have any EFW HX data
% and check for NM/BM1
count_skip = 0;
for cl_id=sc_list
  st_tmp = st;
  tm = []; tm_prev = [];
  
  while st_tmp<st+dt
    irf_log('proc',['Requesting C' num2str(cl_id)...
      ' : ' epoch2iso(st_tmp,1)])
    [t,data] = caa_is_get(DB_S,st_tmp,REQ_INT,cl_id,'efw','FDM'); %#ok<ASGLU>
    if ~isempty(data), tm_cur = data(5,1);
    else
      irf_log('dsrc',['No FDM for C' num2str(cl_id) ...
        ' at ' epoch2iso(st_tmp,1)])
      tm_cur = -1;
    end
    
    % First frame is directly good, probably we started in the middle of
    % good interval
    if isempty(tm_prev) && tm_cur>=0
      tm(end+1,:) = [st_tmp tm_cur]; %#ok<AGROW>
    end
    
    if ~isempty(tm_prev) && tm_prev~=tm_cur && tm_cur>=0
      % We skip MAX_SKIP frames of NM because it is folliwing
      % BM1 or data gap and usually contains junk
      if tm_cur==0 && count_skip<MAX_SKIP
        tm_cur = -2;
        count_skip = count_skip + 1;
        irf_log('proc',['Skipping one NM frame for C' num2str(cl_id)...
          ' at ' epoch2iso(st_tmp,1)])
      else
        tm(end+1,:) = [st_tmp tm_cur]; %#ok<AGROW>
        count_skip = 0;
      end
    else, count_skip = 0;
    end
    
    st_tmp = st_tmp + REQ_INT;
    tm_prev = tm_cur;
  end
  
  % Throw away all modes>1
  if ~isempty(tm)
    ii_out = find(tm(:,2)>1);
    if ~isempty(ii_out)
      irf_log('proc',sprintf('Removing %d intervals TM>1 for C%d',...
        length(ii_out),cl_id))
      tm(ii_out,:) = []; %#ok<AGROW>
    end
  end
  c_eval('tm? = tm;',cl_id)
end
clear tm tm_cur

% Make SC_LIST from SC for which TM is not empty
sc_list = [];
c_eval('if exist(''tm?'',''var''),if ~isempty(tm?), sc_list = [sc_list ?]; end, end')

if isempty(sc_list), irf_log('dsrc','No data'), return, end

for cl_id=sc_list
  cdir = [sdir '/C' num2str(cl_id)];
  if ~exist(cdir, 'dir')
    [SUCCESS,MESSAGE] = mkdir(cdir);
    if SUCCESS, irf_log('save',['Created storage directory ' cdir])
    else, error(MESSAGE)
    end
  end
  c_eval('tm=tm?;',cl_id);
  
  % Split long intervals into SPLIT_INT chunks
  j = 1;
  while 1
    st_tmp = tm(j,1);
    tm_cur = tm(j,2);
    if j==size(tm,1), dt_tmp = st + dt - st_tmp;
    else, dt_tmp = tm(j+1,1) - st_tmp;
    end
    % For BM1 we take SPLIT_INT/3 intervals
    if dt_tmp > SPLIT_INT*4/3*(1-tm_cur*2/3)
      if j==size(tm,1)
        tm(j+1,:) = [tm(j,1)+SPLIT_INT*(1-tm_cur*2/3) tm(j,2)]; %#ok<AGROW>
      else
        tm(j+1:end+1,:) = [tm(j,1)+SPLIT_INT*(1-tm_cur*2/3) tm(j,2); tm(j+1:end,:)];
      end
    end
    if j>=size(tm,1), break
    else, j = j + 1;
    end
  end
  
  % Read NS_OPS database
  clear ns_ops
  ns_ops = c_ctl('get',cl_id,'ns_ops');
  if isempty(ns_ops)
    c_ctl('load_ns_ops',[DP_S '/caa-control'])
    ns_ops = c_ctl('get',cl_id,'ns_ops');
  end
  if isempty(ns_ops), error(['cannot get NS_OPS for C' num2str(cl_id)]), end
  
  % Split intervals if a non-blanking nsops starts in the interval
  ii = find( ns_ops(:,1)>st & ns_ops(:,1)<st+dt);
  problem_list=[caa_str2errid('hxonly') caa_str2errid('bad_bias') caa_str2errid('bad_hx') caa_str2errid('bad_lx')];
  c_eval('problem_list=[problem_list caa_str2errid(''no_p?'')];',[1 2 3 4 12 32 34]);
  if size(ii)>0
    for j=1:size(ii)
      if any(ns_ops(ii(j),4)==problem_list)
        irf_log('proc',['NSOPS interval start/stop. Splitting interval for C' num2str(cl_id)...
          ' at ' epoch2iso(ns_ops(ii(j),1))])
        tm=sort([tm' [ns_ops(ii(j),1) 0]'],2)';
      end
    end
  end
  % Split intervals if a non-blanking nsops stops in the interval
  ii = find(ns_ops(:,1)+ns_ops(:,2)>st & ns_ops(:,1)+ns_ops(:,2)<st+dt);
  if size(ii)>0
    for j=1:size(ii)
      if any(ns_ops(ii(j),4)==problem_list)
        irf_log('proc',['NSOPS interval start/stop. Splitting interval for C' num2str(cl_id)...
          ' at ' epoch2iso(ns_ops(ii(j),1))])
        tm=sort([tm' [ns_ops(ii(j),1)+ns_ops(ii(j),2) 0]'],2)';
      end
    end
  end
  
  for inter=1:size(tm,1)
    t1 = tm(inter,1);
    if inter==size(tm,1), dt1 = st +dt -t1;
    else, dt1 = tm(inter+1,1) -t1;
    end
    
    % Disregard bad NS_OPS intervals directly here
    [st_nsops, dt_nsops] = caa_ns_ops_int(t1,dt1,ns_ops); %#ok<NASGU>
    if isempty(st_nsops)
      irf_log('proc','Bad interval according to NS_OPS')
      irf_log('proc',['C' num2str(cl_id) ' skipping ' ...
        epoch2iso(t1,1) ' -- ' epoch2iso(t1+dt1,1)])
      continue
    end
    
    % Intervals shorter then 300 sec and which are not at the
    % beginning of the entire interval are considered bad,
    % as they are usually signatures of hacked data
    if inter~=1 && inter~=size(tm,1) && dt1<300
      irf_log('proc',['C' num2str(cl_id) ' skipping short int ' ...
        epoch2iso(t1,1) ' -- ' epoch2iso(t1+dt1,1)])
      continue
    else
      irf_log('proc',['C' num2str(cl_id) ' interval ' ...
        epoch2iso(t1,1) ' -- ' epoch2iso(t1+dt1,1)])
    end
    
    % Get the data
    c_get_batch(t1,dt1,'db',DB_S,'sc_list',cl_id,'sdir',cdir,'vars',srcvars,'check_caa_sh_interval',1,'noproc')
  end
end

