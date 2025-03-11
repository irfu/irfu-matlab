function out=c_ctl(varargin)
%C_CTL Cluster control
%
% c_ctl('init',[dir])
% c_ctl('get',cl_id,ctl_name)
% c_ctl(cl_id,ctl_name) % equvalent to get
% c_ctl('list',[sc_list])
% c_ctl('load_ns_ops',[dir])
% c_ctl('load_man_int',[dir])
% c_ctl('save',[dir])
% c_ctl('set',[sc_list],ctl_name,ctl_val)
% c_ctl(sc_list,ctl_name,value,[ctl_name1,value1...]) % equvalent to set
% c_ctl('load_hk_cal')
% c_ctl('load_aspoc_active')
% c_ctl('load_bad_ib');
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


default_mcctl_path = '.';

if nargin<1, c_ctl_usage, return, end
args = varargin;

if ischar(args{1})
  if strcmp(args{1},'init')
    if nargin>=2
      d = args{2};
      if exist(d,'dir')
        if exist([d '/mcctl.mat'],'file')
          clear global c_ct; global c_ct
          eval(['load ' d '/mcctl.mat'])
        else
          error(['No mcctl.mat in ' d])
        end
      elseif exist(d,'file')
        clear global c_ct; global c_ct
        eval(['load ' d])
      else
        error(['Directory or file ' d ' does not exist'])
      end
    else
      clear global c_ct; global c_ct
      if exist([default_mcctl_path '/mcctl.mat'],'file')
        eval(['load ' default_mcctl_path '/mcctl.mat'])
      else
        % init
        % NOTE: this is the place to initialize defaults
        def_ct.ns_ops = [];
        def_ct.man_int = [];
        def_ct.ang_lim = 10;	% algle limit for E.B=0
        def_ct.rm_whip = 1;		% remove times with WHI pulses
        def_ct.probe_p = 34;	% default probe pair to use
        def_ct.deltaof_max = 1.5;
        % maximum reasonable value of deltaof
        def_ct.deltaof_sdev_max = 2;
        % delta offsets we remove points which
        % are > deltaof_sdev_max*sdev

        % DSI offsets are moved to c_efw_dsi_off
        % use c_ctl only if you want to override the default value
        def_ct.dsiof = [];

        def_ct.ibias = zeros(256,5,'single');
        def_ct.puck  = zeros(256,5,'single');
        def_ct.guard = zeros(256,5,'single');

        def_ct.aspoc = NaN;

        def_ct.badib = NaN;

        c_ct{1} = def_ct;
        c_ct{2} = def_ct;
        c_ct{3} = def_ct;
        c_ct{4} = def_ct;
        clear def_ct

        % cell number 5 has global settings
        % this cell must be accessed as SC # 0
        if ismac || isunix
          def_ct.isdat_db = 'db.irfu.se:0';
        else
          def_ct.isdat_db = 'NO ISDAT ON WINDOWS';
        end
        if ismac
          def_ct.data_path = '/Volumes/cluster';
        elseif isunix
          def_ct.data_path = '/data/cluster';
        else
          def_ct.data_path = '';
        end
        def_ct.caa_mode = 0;
        c_ct{5} = def_ct;
      end
    end

  elseif strcmp(args{1},'get')
    if nargin<3, error('get: must be c_ctl(''get'',cl_id,''ct_name''))'), end
    cl_id = args{2}; if cl_id==0, cl_id = 5; end
    c = args{3};

    global c_ct
    if isempty(c_ct)
      irf_log('fcal','CTL is not initialized. Initializing...')
      c_ctl('init')
      global c_ct
    end

    if isfield(c_ct{cl_id},c)
      if nargout>0
        eval(['out=c_ct{cl_id}.' c ';'])
      else
        if cl_id < 5, disp(['C' num2str(cl_id) '->' c ':'])
        else, disp(['GLOBAL->' c ':'])
        end
        eval(['disp(c_ct{cl_id}.' c ');'])
      end
    else
      if nargout>0, out=[]; end
      irf_log('fcal',['unknown ctl: ' c])
    end
  elseif strcmp(args{1},'load_hk_cal')

    global c_ct
    if isempty(c_ct)
      irf_log('fcal','CTL is not initialized. Initializing...')
      c_ctl('init')
      global c_ct
    end

    % Read from (/Volumes)/data/cluster/cal

    [c_ct{1}.ibias, c_ct{1}.puck, c_ct{1}.guard] = readhkcalmatrix('C1_CT_EFW_20001128_V002.cal');
    [c_ct{2}.ibias, c_ct{2}.puck, c_ct{2}.guard] = readhkcalmatrix('C2_CT_EFW_20001128_V002.cal');
    [c_ct{3}.ibias, c_ct{3}.puck, c_ct{3}.guard] = readhkcalmatrix('C3_CT_EFW_20001128_V002.cal');
    [c_ct{4}.ibias, c_ct{4}.puck, c_ct{4}.guard] = readhkcalmatrix('C4_CT_EFW_20001128_V002.cal');
  elseif strcmp(args{1},'load_aspoc_active')

    global c_ct
    if isempty(c_ct)
      irf_log('fcal','CTL is not initialized. Initializing...')
      c_ctl('init')
      global c_ct
    end

    if nargin>1, d = args{2};
    else, d = '.';
    end

    c_ct{1}.aspoc = readaspocactive('/C1_CP_ASP_ACTIVE__20010101_000000_20100101_000000_V081030.cef', d);
    c_ct{2}.aspoc = readaspocactive('/C2_CP_ASP_ACTIVE__20010101_000000_20100101_000000_V081030.cef', d);
    c_ct{3}.aspoc = readaspocactive('/C3_CP_ASP_ACTIVE__20010101_000000_20100101_000000_V081030.cef', d);
    c_ct{4}.aspoc = readaspocactive('/C4_CP_ASP_ACTIVE__20010101_000000_20100101_000000_V081030.cef', d);
  elseif strcmp(args{1},'load_bad_ib')

    global c_ct
    if isempty(c_ct)
      irf_log('fcal','CTL is not initialized. Initializing...')
      c_ctl('init')
      global c_ct
    end

    if nargin>1, d = args{2};
    else, d = '.';
    end

    badiblist=ibfn2epoch('/BAD_IB_L2_LIST_C1-4.txt', d);
    for i=1:4
      c_ct{i}.badib=badiblist{i};
      %size(c_ct{i}.badib)
    end
  elseif strcmp(args{1},'load_ns_ops')
    global c_ct
    if isempty(c_ct)
      irf_log('fcal','CTL is not initialized. Initializing...')
      c_ctl('init')
      global c_ct
    end

    if nargin>1, d = args{2};
    else, d = '.';
    end

    for j=1:4
      try
        f_name = [d '/ns_ops_c' num2str(j) '.dat'];
        if exist(f_name,'file')
          c_ct{j}.ns_ops = load(f_name,'-ascii');

          % remove lines with undefined dt
          c_ct{j}.ns_ops(find(c_ct{j}.ns_ops(:,2)==-157),:) = [];
        else, irf_log('load',['file ' f_name ' not found'])
        end
      catch
        disp(lasterr)
      end
    end

  elseif strcmp(args{1},'load_man_int')
    global c_ct
    if isempty(c_ct)
      irf_log('fcal','CTL is not initialized. Initializing...')
      c_ctl('init')
      global c_ct
    end

    if nargin>1, d = args{2};
    else, d = '.';
    end

    for j=1:4
      try
        f_name = [d '/QRecord_c' num2str(j) '.dat'];
        if exist(f_name,'file')
          c_ct{j}.man_int = load(f_name,'-ascii');

          % remove lines with undefined dt
          c_ct{j}.man_int(find(c_ct{j}.man_int(:,2)==-157),:) = [];      % TODO: Is this valid for man_int as well?
        else, irf_log('load',['file ' f_name ' not found'])
        end
      catch
        disp(lasterr)
      end
    end

  elseif strcmp(args{1},'list')
    if nargin>=2
      sc_list = args{2};
      ii = find(sc_list==0);
      if ~isempty(ii), sc_list(ii) = 5; end
    else, sc_list=1:5;
    end
    global c_ct
    if isempty(c_ct), disp('CTL is not initialized.'), return, end
    for cl_id=sc_list
      c = fieldnames(c_ct{cl_id});
      if ~isempty(c)
        for j=1:length(c)
          if cl_id < 5, disp(['C' num2str(cl_id) '->' c{j} ':'])
          else, disp(['GLOBAL->' c{j} ':'])
          end
          eval(['disp(c_ct{cl_id}.' c{j} ');'])
        end
      end
    end

  elseif strcmp(args{1},'save')
    global c_ct
    if isempty(c_ct), disp('CTL is not initialized.'), return, end

    if nargin>1
      d = args{2};
      if exist(d,'dir')
        disp(['Saving ' d '/mcctl.mat'])
        eval(['save -MAT ' d '/mcctl.mat c_ct'])
      else
        disp(['Saving ' d])
        eval(['save -MAT ' d ' c_ct'])
      end
    else
      disp('Saving mcctl.mat')
      save -MAT mcctl.mat c_ct
    end

  elseif strcmp(args{1},'set')
    if nargin<3, error('set: must be c_ctl(''set'',''ctl_name'',value))'), end
    if isnumeric(args{2})
      sc_list = args{2};
      ii = find(sc_list==0);
      if ~isempty(ii), sc_list(ii) = 5; end
      c = args{3};
      c_val = args{4};
    else
      sc_list = 1:4;
      c = args{2};
      c_val = args{3};
    end
    if ~ischar(c), error('ctl_name must be a string'), end
    global c_ct
    if isempty(c_ct), disp('CTL is not initialized.'), return, end
    for cl_id=sc_list
      try
        eval(['c_ct{cl_id}.' c '=c_val;'])
      catch
        disp(lasterr)
        error('bad option')
      end
      if cl_id < 5, disp(['C' num2str(cl_id) '->' c ':'])
      else, disp(['GLOBAL->' c ':'])
      end
      eval(['disp(c_ct{cl_id}.' c ');'])
    end

  else
    error('Invalid argument')
  end
elseif isnumeric(args{1})
  sc_list = args{1};
  ii = find(sc_list==5);
  if ~isempty(ii), sc_list(ii) = 5; end

  if nargin>2, have_options = 1; args = args(2:end);
  elseif nargin==2
    have_options = 0;
    if nargout>0, out=c_ctl('get',sc_list,args{2});
    else, c_ctl('get',sc_list,args{2});
    end
  else, have_options = 0;
  end

  while have_options
    if length(args)>1
      if ischar(args{1})
        c_ctl('set',sc_list,args{1},args{2})
      else
        error('option must be a string')
      end
      if length(args) >= 2
        args = args(3:end);
        if isempty(args), break, end
      else, break
      end
    else
      disp('Usage: c_ctl(sc_list,''ctl_name'',value)')
      break
    end
  end
else
  error('Invalid argument')
end

function ret = findhkcalmatrix( fid, searchstr )
% Find string in hk calib matrix
ret = -2;
tline = fgetl(fid);
while ischar(tline)
  %        disp(tline)
  if strncmpi(tline,searchstr,length(searchstr))
    ret = -1;
    break;
  end
  tline = fgetl(fid);
end
while ischar(tline)
  %        disp(tline)
  if strncmpi(tline,'# step',6) % find first data line
    ret = 0;
    break;
  end
  tline = fgetl(fid);
end

function [ibias, puck, guard] = readhkcalmatrix( filen )
% Read hk calib matrix from file

if ismac
  datapath = '/Volumes/cluster/cal/';
else
  datapath = '/data/cluster/cal/';
end

fid = fopen([ datapath filen ], 'r');

if fid >= 0
  ibias=zeros(256,5,'single');
  ret = findhkcalmatrix(fid, 'QTY          IBIAS1');
  if ret ~= 0
    irf_log('load',['BIAS hk cal matrix not found in ' datapath filen]);
  else
    ib = textscan(fid,' %d16 %f32 %f32 %f32 %f32',255);
    for i=1:5
      ibias(2:end,i)=ib{i};
      if i > 1 % 1st value computed as in isdat ./server/Wec/Efw/calib_read.c
        ibias(1,i)=2*ibias(2,i)-ibias(3,i);
      end
    end
  end

  puck=zeros(256,5,'single');
  ret = findhkcalmatrix(fid, 'QTY          PUCK1');
  if ret ~= 0
    irf_log('load',['PUCK hk cal matrix not found in ' datapath filen]);
  else
    ib = textscan(fid,' %d16 %f32 %f32 %f32 %f32',255);
    for i=1:5
      puck(2:end,i)=ib{i};
      if i > 1 % 1st value computed as in isdat ./server/Wec/Efw/calib_read.c
        puck(1,i)=2*puck(2,i)-puck(3,i);
      end
    end
  end

  guard=zeros(256,5,'single');
  ret = findhkcalmatrix(fid, 'QTY          GUARD1');
  if ret ~= 0
    irf_log('load',['GUARD hk cal matrix not found in ' datapath filen]);
  else
    ib = textscan(fid,' %d16 %f32 %f32 %f32 %f32',255);
    for i=1:5
      guard(2:end,i)=ib{i};
      if i > 1 % 1st value computed as in isdat ./server/Wec/Efw/calib_read.c
        guard(1,i)=2*guard(2,i)-guard(3,i);
      end
    end
  end

  fclose(fid);
else
  irf_log('load',['File ' datapath filen ' not found']);
end

function aspa = readaspocactive( filen, datapath )
%READASPOCACTIVE read cef ascii files
%   aa = read_cef_active(file) caa cef file and returns the data.
%   aa(x,1) is the ON date in epoch format
%   aa(x,2) is the OFF date in epoch format
%
%    if ismac
%        datapath = '/Volumes/caa/cef/ASPOC/';
%    else
%        datapath = '/data/caa/cef/ASPOC/';
%    end

fid = fopen([ datapath filen ], 'r');

aspa=[];
cols=0;
if fid==-1
  irf_log('load',['File ' datapath filen ' not found']);
  return
end
s = '';
while not(strncmp(s,'DATA_UNTIL',10))
  s = fgets(fid);
end

fposmem=ftell(fid);
s = fgets(fid);
if s(1) == '!' % check for no data
  fclose(fid);
  %    disp(['Reading: ' file ' done. No data.']);
  return;
end
% Count columns
cols=1;
for i=1:size(s,2)
  if s(i)=='/'
    cols=cols+1;
  end
end

fseek(fid,fposmem,'bof');

if cols==2
  C = textscan(fid, '%24s/%24s $');
else
  disp('PANIC: Can only read cef files with on/off data.');
  return;
end
fclose(fid);
C1=C{1}(1:end-1);
C2=C{2}(1:end);
for i=1:size(C1,1)
  aspa(i,1) = iso2epoch(C1{i});
  aspa(i,2) = iso2epoch(C2{i});
end

function timeclid = ibfn2epoch( filen, datapath )
% Make an epoch list of burst start time and cluster ID from
% a text file with burst file names from /(data|Volumes)/cluster/burst.

%    if ismac
%        datapath = '/Volumes/caa/cef/iburst/';
%    else
%        datapath = '/data/caa/cef/iburst/';
%    end

fid = fopen([ datapath filen ],'r'); %Open the text file that contain the information about the burst.
if fid==-1
  irf_log('load',['File ' datapath filen ' not found']);
  return
end
i=0;
cnt=zeros(1,4);
timeclid={ [] [] [] [] };
s='Empty';
while ~feof(fid)
  tline = fgetl(fid);
  if length(tline)>=17
    s = tline(1:17);
    cl = str2double(s(end));
    if cl<1 || cl>4
      error(['Spacecraft ID must be 1..4. Value: ' num2str(cl) ]);
    end
    full_time = iso2epoch(['20' s(1:2) '-' s(3:4) '-' s(5:6) 'T' s(7:8) ':' s(9:10) ':' s(11:12) 'Z']);
    cnt(cl)=cnt(cl)+1;
    timeclid{cl}(cnt(cl),1)=full_time;
  end
end
%cnt
irf_log('proc',['Last bad L2 iburst file name read: ' s]);
fclose(fid);


function c_ctl_usage
disp('Usage:')
disp('  c_ctl(''init'')')
disp('  c_ctl(''init'',''/path/to/mcctl.mat/'')')
disp('  c_ctl(''init'',''/path/to/alternative_mcctl.mat'')')
disp('  c_ctl(''load_ns_ops'')')
disp('  c_ctl(''load_ns_ops'',''/path/to/ns_ops_cN.dat/'')')
disp('  c_ctl(''load_man_int'')')
disp('  c_ctl(''load_man_int'',''/path/to/QRecord_cN.dat/'')')
disp('  c_ctl(''save'')')
disp('  c_ctl(''save'',''/path/to/mcctl.mat/'')')
disp('  c_ctl(''save'',''/path/to/alternative_mcctl.mat'')')
disp('  c_ctl(''sc_list'',''ctl'',value)')
disp('  c_ctl(''load_hk_cal'')')
disp('  c_ctl(''load_aspoc_active'')')
disp('  c_ctl(''load_bad_ib'')')
