function c_get_batch(st,dt,varargin)
%C_GET_BATCH prepare data from all SC: do all steps until spinfits and despin
%
% c_get_batch(start_time,dt/stop_time,[sc_list],[options ...])
% Input:
% start_time - ISDAT epoch/ISO time string
% dt/stop_time - length of time interval in sec
% sc_list - list of SC [optional]
% Options: go in pair 'option', value
% 'sp' - storage directory;
%   // default: '.'
% 'sdir' - main storage directory, storage directory is constructed from SDIR
%   and start_time: SDIR/YYYYMMDD_hhmm;
% 'dp' - storage directory;
%   // default: '/data/cluster'
% 'db' - ISDAT database;
%   // default: from c_ctl
% 'sc_list' - list of SC;
%   // default: 1:4
% 'vars' - variables to get data for (see help ClusterDB/getData)
%   supplied as a string separated by '|' or as a cell array;
%   // default: {'tmode','fdm','efwt','ibias','p','e','a','sax',...
%   'r','v','b','edi','ncis','vcis','bfgm'}
%   'ncis','vcis','vce','bfgm'}
%   // + {'dies','die','brs','br'} which are always added. Use 'noproc'
%   to skip them.
% 'varsproc' - variables to get data for (see help ClusterProc/getData)
% 'nosrc'  - do not run ClusterDB/getData;
% 'noproc' - do not run ClusterProc/getData for {'dies','die','brs','br'};
% 'extrav' - extra variables in addition to default;
% 'swmode'  - do solar wind wake cleaning for E-field
% 'cdb' - ClusterDB object;
% ++ extra options to be passwed to ClusterProc/getData
%
% Examples:
% c_get_batch(toepoch([2002 03 04 10 00 00]),30*60,...
% 'sp','/home/yuri/caa-data/20020304')
%   % load all data to /home/yuri/caa-data/20020304
% c_get_batch(toepoch([2002 03 04 10 00 00]),30*60,...
% 'sp','/home/yuri/caa-data/20020304','vars',{'e','p'})
%   % load only 'e' and 'p'
% c_get_batch(toepoch([2002 03 04 10 00 00]),30*60,...
% 'sp','/home/yuri/caa-data/20020304','vars','e','nosrc','withwhip')
%   % to recompute E and pass 'withwhip' to ClusterProc/getData
%
% See also ClusterDB/getData, ClusterProc/getData
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

persistent st_vector dt_sec;

narginchk(0 ,15)

sc_list = 1:4;

if nargin>2, have_options = 1; args = varargin;
elseif nargin==0
  have_options=0;
  st_vector=irf_ask('Start time [%]>','st_vector',[2002 01 01 0 0 0]);
  st=toepoch(st_vector);
  dt_sec=irf_ask('Time interval in seconds [%]>','dt_sec',60);
  dt=dt_sec;
elseif nargin==1
  error('do not udnerstand what to do with the argument, see help');
else, have_options = 0;
end

sp = '.';
db = c_ctl(0,'isdat_db');
dp = c_ctl(0,'data_path');
cdb = '';
vars = '';
varsXtra = '';
varsProc = '';
argsProc = '';
dosrc = 1;
doproc = 1;
sw_mode = 0;
sdir_st = '';



if have_options
  if isnumeric(args{1})
    sc_list = args{1};
    if length(args)>1, args = args(2:end);
    else, have_options = 0;
    end
  end
end

while have_options
  l = 2;
  if length(args)>=1
    switch(args{1})
      case 'sp'
        if ischar(args{2}), sp = args{2};
        else, irf_log('fcal','SP must be string')
        end
      case 'sdir'
        if ischar(args{2}), sdir_st = args{2};
        else, irf_log('fcal','SDIR must be string')
        end
      case 'dp'
        if ischar(args{2}), dp = args{2};
        else, irf_log('fcal','DP must be string')
        end
      case 'db'
        if ischar(args{2}), db = args{2};
        else, irf_log('fcal','DB must be string')
        end
      case 'sc_list'
        if isnumeric(args{2}), sc_list = args{2};
        else, irf_log('fcal','SC_LIST must be numeric')
        end
      case 'swmode'
        sw_mode = 1;
      case 'vars'
        if ischar(args{2})
          vars = {};
          p = tokenize(args{2},'|');
          for i=1:length(p), vars(length(vars)+1) = p(i); end %#ok<AGROW>
        elseif iscell(args{2}), vars = args{2};
        else
          irf_log('fcal','VARS must be either string or cell array')
        end
      case 'varsproc'
        if ischar(args{2})
          varsProc = {};
          p = tokenize(args{2},'|');
          for i=1:length(p), varsProc(length(varsProc)+1) = p(i); end %#ok<AGROW>
        elseif iscell(args{2}), varsProc = args{2};
        else
          irf_log('fcal',...
            'VARSPROC must be either string or cell array')
        end
      case 'extrav'
        if ischar(args{2})
          varsXtra = tokenize(args{2},'|');
        elseif iscell(args{2}), varsXtra = args{2};
        else
          irf_log('fcal',...
            'EXTRAV must be either string or cell array')
        end
      case 'cdb'
        if (isa(args{2},'ClusterDB')), cdb = args{2};
        else
          irf_log('fcal','CDB must be a ClusterDB object')
        end
      case 'nosrc'
        dosrc = 0; l = 1;
      case 'noproc'
        doproc = 0; l = 1;
      case 'check_caa_sh_interval'
        argsProc = [{'check_caa_sh_interval'} argsProc]; %#ok<AGROW>
      case 'ec_extraparams'
        argsProc = [argsProc {'ec_extraparams'} args(2)]; %#ok<AGROW>
      otherwise
        irf_log('fcal',...
          ['Option ''' args{1}...
          ''' not recognized. Pass the rest to getData'])
        argsProc = args;
        break
    end
    if length(args) > l, args = args(l+1:end);
    else, break
    end
  else
    error('caa:wrongArgType','use c_get_batch(..,''option'',''value'')')
  end
end

if dosrc==1, [st,dt] = irf_stdt(st,dt); end
if ~isempty(sdir_st), sp = [sdir_st '/' irf_fname(st)]; end

if isempty(cdb), cdb = ClusterDB(db,dp,sp); end

if isempty(vars) && isempty(varsProc)
  vars = {'tmode','fdm','efwt','ibias','p','e','a','sax',...
    'r','v','b','edi','ncis','vcis','bfgm'};
end

if ~isempty(varsXtra), vars = [vars varsXtra]; end

if ~isempty(vars)
  if L_find(vars,{'e','p'})
    if isempty(L_find(vars,'ibias')), vars = [{'ibias'} vars]; end
    if isempty(L_find(vars,'efwt')), vars = [{'efwt'} vars]; end
    if isempty(L_find(vars,'fdm')), vars = [{'fdm'} vars]; end
  end
  if ~isempty(L_find(vars,{'e','bsc'})) &&  isempty(L_find(vars,'tmode'))
    vars = [{'tmode'} vars];
  end

  if dosrc
    for cl_id=sc_list
      for k=1:length(vars), getData(cdb,st,dt,cl_id,vars{k}); end
    end
  end

  if doproc && isempty(varsProc)
    if L_find(vars,{'e','p'})
      varsProc = [{'whip','sweep','bdump'} varsProc];
    end
    if L_find(vars,{'ibias','efwt'}), varsProc = [varsProc {'badbias'}]; end
    if L_find(vars,'p'), varsProc = [varsProc {'probesa','p','ps'}]; end
    if L_find(vars,'e')
      if sw_mode
        varsProc = [varsProc {'ec','dies','die'}];
        argsProc = [{'correct_sw_wake','rmwhip'} argsProc];
      else
        varsProc = [varsProc {'dies','die'}];
      end
    end
    if L_find(vars,{'b','bfgm'}), varsProc = [varsProc {'brs','br'}]; end
    if L_find(vars,{'bsc'}), varsProc = [varsProc {'dibsc'}]; end
    if L_find(vars,'edi'), varsProc = [varsProc {'edi'}]; end
    if L_find(vars,{'pburst','eburst'}), varsProc = [varsProc {'dieburst'}]; end
    if L_find(vars,'vcis'), varsProc = [varsProc {'vce'}]; end
  end
end

if ~isempty(varsProc) && doproc
  cp=ClusterProc(get(cdb,'sp'));
  for cl_id=sc_list
    for k=1:length(varsProc)
      proc_options = [{cp} {cl_id} varsProc(k) argsProc];
      getData(proc_options{:});
    end
  end
end

function ii = L_find(list,s_list)
ii = [];
if isempty(list), return, end
if ischar(s_list)
  % fast search
  for j=1:length(list)
    if strcmp(list{j},s_list), ii = j; return, end
  end
else
  for k=1:length(s_list)
    for j=1:length(list)
      if strcmp(list{j},s_list{k}), ii = [ii j]; break, end %#ok<AGROW>
    end
  end
end
