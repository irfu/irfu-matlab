function ret = caa_get_bursts(filename, plot_flag)
%caa_get_bursts(filename) produce and plot Cluster burst data from the raw data
%
% data = caa_get_bursts(filename, plot_flag)
%
% Input:
%   filename   : burst data file name
%   plot_flag  : generate plot 0=off(default) 1=on 2=on+local data (not CAA)
%
% mEFWburstTM.mat, mEFWburstR.mat & mEFWburst.mat are saved in
% data/caa/l1 working directory.
% Burst files are read from directory /data/cluster/burst/
% Plot files are saved in directory $HOME/figures/
%
% Example:
%       caa_get_bursts('070831101839we.03')
%

narginchk(1,2);

flag_local = 0;
if nargin < 2
  plot_flag = 0;
elseif plot_flag==2
  flag_local = 1;
end

plot_save=1;
if flag_local
  plotpath='./';
else
  plotpath=[getenv('HOME') '/figures/'];
end
flag_save = 1;

DP = c_ctl(0,'data_path');
DB = c_ctl(0,'isdat_db');
B_DT = 300;
B_DELTA = 60;

cl_id=str2double(filename(end)); %Get the satellite number
if ( cl_id < 1 || cl_id > 4 )
  error(['Wrong spacecraft ID read from filename. ' filename ])
end
fname=irf_ssub([plotpath 'p?-c!'],filename(1:12),cl_id); %Sets the name that will be used to save the plots
fnshort=filename;
s=filename;
full_time = iso2epoch(['20' s(1:2) '-' s(3:4) '-' s(5:6) 'T' s(7:8) ':' s(9:10) ':' s(11:12) 'Z']);
start_time=full_time;
st=full_time;

if flag_local
  old_pwd = '.';
  [loc_st,loc_dt] = caa_read_interval;
  if isempty(loc_st)
    irf_log('proc','No Cluster data in the current directory (.interval not found)');
    ret=2;
    return;
  end
  loc_st = iso2epoch(loc_st);
  if start_time < loc_st || loc_st > start_time + loc_st
    irf_log('proc',['Current directory(' irf_disp_iso_range(loc_st+[0 loc_dt],1) ')'])
    irf_log('proc',['does not match the burts start time (' epoch2iso(start_time,1) ')']);
    ret=2;
    return;
  end
else
  dirs = caa_get_subdirs(st, 90, cl_id);
  if isempty(dirs)
    irf_log('proc',['Can not find L1 data dir for ' s]);
    ret=1;
    return;
  end
  found=false;
  for i=size(dirs,2):-1:1 % find start time directory
    d=dirs{i}(end-12:end);
    dtime=iso2epoch([d(1:4) '-' d(5:6) '-' d(7:8) 'T' d(10:11) ':' d(12:13) ':00Z']);
    if dtime<=start_time
      found=true;
      break;
    end
  end
  if ~found
    irf_log('proc','iburst start time does not match any L1 data dir');
    ret=2;
    return;
  end
  old_pwd = pwd;
  cd(dirs{i})
end

varsb = c_efw_burst_param([DP '/burst/' filename]);
varsbsize = length(varsb);

if 0    % Exit if no V43M/H or V12H parameter
  found=false; %#ok<UNRCH>
  for i=1:varsbsize
    if strcmp(varsb(i),'V43M') || strcmp(varsb(i),'V43H') || strcmp(varsb(i),'V12H')
      found=true;
      break;
    end
  end
  if ~found
    irf_log('proc','iburst special no V43M/H or V12H continue');
    ret=-1;
    cd(old_pwd);
    return;
  end
end

if ~flag_local
  %Remove old files
  fn={'mEFWburstTM.mat' 'mEFWburstR.mat' 'mEFWburst.mat'};
  for i=1:size(fn,2)
    if exist(fn{i},'file')
      delete(fn{i});
    end
  end
end

for out = 1:varsbsize
  [field,sen,filter] = get_ib_props(varsb{out});
  instrument = 'efw';

  [t,data] = caa_is_get(DB,st-B_DELTA,B_DT,cl_id,instrument,field,...
    sen,filter,'burst','tm');

  if (out==1) % Create data matrix for t and all 8 possible variables
    if size(t,1)<3 || size(data,1)<3 % sanity check
      irf_log('proc','No usable burst data');
      ret=3;
      return;
    end
    data8 = NaN(size(data,1),9);

    start_satt = c_efw_burst_chkt(DB,[DP '/burst/' filename]);
    if isempty(start_satt)
      irf_log('dsrc','burst start time was not corrected')
    else
      err_t = t(1) - start_satt;
      if abs(err_t) > 5
        irf_log('dsrc',['burst start offset ' ...
          num2str(err_t) ' sec is too big!'])
        irf_log('dsrc','burst start time was not corrected')
      else
        irf_log('dsrc',['burst start time was corrected by ' ...
          num2str(err_t) ' sec'])
        t = t - err_t;
      end
    end
    data8(:,1)=t;   % corrected time
  end
  if isempty(data)
    irf_log('proc',['Bad burst data column: ' num2str(out)]);
    ret=4;
    return;
  else
    data8(:,out+1)=data;
  end
end
clear data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Checking the order of the data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if varsbsize>=2
  can_check_order = 1;
  [ok,pha] = c_load('Atwo?',cl_id);
  if ~ok || isempty(pha)
    irf_log('load',...
      irf_ssub('No/empty Atwo? in mA. Use getData(CDB,...,cl_id,''a'')',cl_id))
    can_check_order = 0;
  end

  if can_check_order
    % See which single-ended probes we got in the burst, e.g. V1H, V2H...
    probes = [];
    for i=1:length(varsb)
      if varsb{i}(1)=='V' && ischar(varsb{i}(3))
        probes = [probes str2double(varsb{i}(2))]; %#ok<AGROW>
      end
    end
    if isempty(probes)
      can_check_order = 0;
      irf_log('proc','Cannot check burst order: no single-ended probes')
    end
  end

  if can_check_order
    % See which NM data we may use
    probep = [];
    if any(probes==3) && any(probes==4), probep = 34; end
    if any(probes==1) && any(probes==2), probep = [probep 12]; end
    if any(probes==3) && any(probes==2), probep = [probep 32]; end
    if isempty(probep)
      can_check_order = 0;
      irf_log('proc','Cannot check burst order: no useful probe pairs')
    end
  end

  if can_check_order
    % Load NM data
    for i=1:length(probep)
      [ok,nmdata] = c_load(irf_ssub('wE?p!',cl_id,probep(i)));
      if ok && ~isempty(nmdata)
        ref_probep = probep(i);
        irf_log('proc',sprintf('Using p%d NM data as reference',probep(i)))
        break
      end
    end
    if isempty(nmdata)
      can_check_order = 0;
      irf_log('proc','Cannot check burst order: no reference NM data')
    end
  end

  if can_check_order
    data8ord = c_efw_burst_order_signals(data8, varsb, nmdata, pha, ref_probep);
  else
    data8ord = data8;
  end
else
  data8ord=data8; % assume no order change
end

% save raw ordered data
save_file = './mEFWburstTM.mat';
save_list='';
burst_info=sprintf('%s ',varsb{:}); %#ok<NASGU>
eval(irf_ssub(['ib?_info=burst_info(1:end-1);' 'save_list=[save_list ''ib?_info ''];'],cl_id));

eval(irf_ssub(['iburst?=data8ord;' 'save_list=[save_list ''iburst? ''];'],cl_id));
if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
  irf_log('save',[save_list ' -> ' save_file])
  if exist(save_file,'file')
    eval(['save -append ' save_file ' ' save_list]);
  else
    eval(['save ' save_file ' ' save_list]);
  end
end
% filter data
% calib data

[data8ordfc, ix] = caa_identify_ib_spikes(data8ord);

BSCpos=zeros(1,3);
BSCcnt=0;
for i=1:varsbsize
  if varsb{i}(1)=='B' % BP data
    % BP12 factor is unkrown. These data do not have an easy physical
    % meaning. typically used as a place holder.
    %data8ordfc(:,i+1)=data8ordfc(:,i+1)*BPFACTOR;
  elseif varsb{i}(1)=='S' % SC data
    xyzord=double(varsb{i}(3))-87; % X=1
    BSCpos(xyzord)=i;
    data8ordfc(:,i+1) = -data8ordfc(:,i+1)/7000;
    c_efw_burst_bsc_tf(data8ordfc(:,[1 BSCpos+1]),cl_id,xyzord);
    BSCcnt=BSCcnt+1;
  elseif varsb{i}(1)=='V'
    if length(varsb{i})>3 % Differential signals are special
      probe = str2double(varsb{i}(2:3));
      if probe==12 || probe==43
        data8ordfc(:,i+1) = data8ordfc(:,i+1)*0.00212/0.088; % Convert to mV/m
      else
        % XXX: do we ever have asymmetric, e.g. p23 ?
        error('Strange probe configuration. Update the code!!')
      end
      filter = varsb{i}(end);
      if filter=='H', filter = 'BP'; end% It is a special case
      data8ordfc(:,[1 i+1]) = c_efw_invert_tf(data8ordfc(:,[1 i+1]),filter);
    else
      data8ordfc(:,i+1) = data8ordfc(:,i+1)*0.00212;
      data8ordfc(:,[1 i+1]) = c_efw_invert_tf(data8ordfc(:,[1 i+1]),...
        varsb{i}(end));
    end
  else
    error(['Unknown ib data type: ' varsb{i}]);
  end
end
data8ordfc(ix,2:end) = NaN; % Set spikes to NaNs

% Save data
save_file = './mEFWburstR.mat';
save_list='';
if BSCcnt % BSC
  % Sanity check. 4 param data with only 2 BSC exists 2011 ex:
  % 110814220936we.03 110816184430we.03 110828095728we.01
  if BSCcnt<3
    irf_log('save','No BSC data save. Less than 3 BSC components.');
  else
    BSCtemp = data8ordfc(:,[1 BSCpos+1]);  %#ok<NASGU>
    c_eval('wBSC4kHz?=BSCtemp;save_list=[save_list '' wBSC4kHz? ''];',cl_id);
    clear BSCtemp
  end
end

for i=1:varsbsize
  if varsb{i}~='V'
    continue
  end
  [~,~,filter,probe] = get_ib_props(varsb{i});

  %%%%%%%%%%%%%%%%%%%%%%%%% PROBE MAGIC %%%%%%%%%%%%%%%%%%%%%%
  start_time = data8ordfc(1,1);
  switch cl_id
    case 1
      if start_time>toepoch([2018 12 10 03 00 16])
        % p3 failure
        if any(probe=='1') || any(probe=='3') || any(probe=='4')
          irf_log('dsrc',sprintf('p1, p3 & p4 are BAD on sc%d',cl_id))
          continue
        end
      elseif start_time>toepoch([2009 10 14 07 00 00]) || ...
          (start_time>toepoch([2009 04 19 00 00 00]) && start_time<toepoch([2009 05 07 00 00 00]))
        % p4 failure
        if any(probe=='1') || any(probe=='4')
          irf_log('dsrc',sprintf('p1 & p4 are BAD on sc%d',cl_id))
          continue
        end
      elseif start_time>toepoch([2001 12 28 03 00 00])
        % p1 failure
        if any(probe=='1')
          irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
          continue
        end
      elseif ( (start_time>=toepoch([2001 04 12 03 00 00]) && start_time<toepoch([2001 04 12 06 00 00])) || ...
          (  start_time>=toepoch([2001 04 14 06 00 00]) && start_time<toepoch([2001 04 16 15 00 00])) || ...
          (  start_time>=toepoch([2001 04 18 03 00 00]) && start_time<toepoch([2001 04 20 09 00 00])) || ...
          (  start_time>=toepoch([2001 04 21 21 00 00]) && start_time<toepoch([2001 04 22 03 00 00])) || ...
          (  start_time>=toepoch([2001 04 23 09 00 00]) && start_time<toepoch([2001 04 23 15 00 00])) )
        % The bias current is a bit too large
        % on p3 and p4 on C1&2 in April 2001.
        % Ignore p3, p4 and p34 and only use p1, p2 and p12.
        % Use only complete 3-hour intervals to keep it simple.
        if any(probe=='3') || any(probe=='4')
          irf_log('dsrc',sprintf('Too high bias current on p3 & p4 sc%d',cl_id));
          continue
        end
      end
    case 2
      if start_time>toepoch([2022 08 23 12 08 0])
        % p3 failure
        if any(probe=='1') || any(probe=='2') || any(probe=='3')
          irf_log('dsrc',sprintf('p1, p2 & p3 is BAD on sc%d',cl_id))
          continue
        end
        elseif start_time>toepoch([2015 10 12 12 00 0])
        % p2 failure
        if any(probe=='1') || any(probe=='2')
          irf_log('dsrc',sprintf('p1 & p2 is BAD on sc%d',cl_id))
          continue
        end
      elseif start_time>=toepoch([2007 05 13 03 23 48])
        % p1 failure
        if any(probe=='1')
          irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
          continue
        end
      end
    case 3
      if start_time>toepoch([2014 11 03 20 58 16.7])
        % p2 failure
        if any(probe=='1') || any(probe=='2') || any(probe=='3')
          irf_log('dsrc',sprintf('p1, p2 & p3 are BAD on sc%d',cl_id));
          continue
        end
      elseif start_time>toepoch([2011 6 01 09 30 0])
        % p3 failure
        if any(probe=='1') || any(probe=='3')
          irf_log('dsrc',sprintf('p1 & p3 are BAD on sc%d',cl_id));
          continue
        end
      elseif start_time>toepoch([2002 07 29 09 06 59 ])
        % p1 failure
        if any(probe=='1')
          irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
          continue
        end
      end
    case 4
      if start_time>=toepoch([2015 02 17 07 30 00]) % 2015-02-17 07:36:30
        % p3 failure
        if any(probe=='3') || any(probe=='4')
          irf_log('dsrc',sprintf('p3 & p4 are BAD on sc%d',cl_id));
          continue
        end
      elseif start_time>=toepoch([2013 07 01 13 30 00]) % 2013 07 01 14 43 44
        % p4 failure
        if any(probe=='4')
          irf_log('dsrc',sprintf('p4 is BAD on sc%d',cl_id));
          continue
        end
      end
  end
  %%%%%%%%%%%%%%%%%%%%%%% END PROBE MAGIC %%%%%%%%%%%%%%%%%%%%

  data=data8ordfc(:,[1 i+1]);
  if length(probe)>1
    % Invert the sign to match the HX data
    data(:,2) = -data(:,2); %#ok<NASGU>
    if strcmp(probe,'43'), probe = '34'; end
    eval(irf_ssub(['wE?!p$=data;' 'save_list=[save_list ''wE?!p$ ''];'],filter,cl_id,probe));
  else
    eval(irf_ssub(['P?!p$=data;' 'save_list=[save_list ''P?!p$ ''];'],filter,cl_id,probe));
  end
end

% Make spike problem time vector
mem=1;
ii=1;
tinter=[];
sz=length(ix);
if sz>1
  for i=2:sz
    if ix(i)-ix(i-1)>1
      tinter(ii,1)=data8ordfc(ix(mem),1); %#ok<AGROW>
      tinter(ii,2)=data8ordfc(ix(i-1),1); %#ok<AGROW>
      ii=ii+1;
      mem=i;
    end
  end
  tinter(ii,1)=data8ordfc(ix(mem),1);
  tinter(ii,2)=data8ordfc(ix(sz),1);
end
if ~isempty(tinter)
  eval(irf_ssub(['SPIKE?=tinter;' 'save_list=[save_list ''SPIKE? ''];'],cl_id));
end
clear tinter ii mem sz

if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
  irf_log('save',[save_list ' -> ' save_file])

  if exist(save_file,'file')
    eval(['save -append ' save_file ' ' save_list]);
  else
    eval(['save ' save_file ' ' save_list]);
  end
end

% make L2
cp = ClusterProc;
if ~flag_local
  getData(cp,cl_id,'whip');
end
getData(cp,cl_id,'pburst');
getData(cp,cl_id,'dieburst');
getData(cp,cl_id,'dibscburst');

if plot_flag
  st=data8ordfc(1,1);
  sp=data8ordfc(end,1);
  if abs(st)==Inf || isnan(st)
    ret=5;
    irf_log('proc','Cannot plot IB. Bad time vector.')
    return;
  end
  if ~flag_local
    getData(ClusterDB(DB,DP,'.'),st-B_DELTA,B_DT,cl_id,'bfgm');
  end

  clf;
  dt2=5;
  st_int=st-dt2;
  st_int2=sp-st+2*dt2;

  summaryPlot(cp,cl_id,'fullb','ib','st',st_int,'dt',st_int2,'vars',char([fnshort varsb]));

  if plot_save
    orient landscape
    irf_print_fig(fname,'pdf')
  end
end

ret=0;
cd(old_pwd);
end

function [field,sen,filter,probe] = get_ib_props(ibvarstr)
% get filter from ib string
field='E';
switch ibvarstr(1)
  case 'B'
    probe=ibvarstr(3:4);
    sen = irf_ssub('p?',probe);
    filter='bp';
  case 'S'
    field='dB';
    probe=lower(ibvarstr(3));
    sen=probe;
    filter='4kHz';
  case 'V'
    if length(ibvarstr) > 3
      probe = ibvarstr(2:3);
      switch probe
        case '12'
          sen = 'p12';
        case '43'
          sen = 'p34'; % for isdat
        otherwise
          error('Bad probe value "%s". Must be 12 or 43',probe)
      end
    else
      probe = ibvarstr(2);
      sen = irf_ssub('p?',probe);
    end
    switch ibvarstr(end)
      case 'U'
        filter='32kHz';
      case 'H'
        if length(probe)==1, filter='4kHz';
        else, filter='8kHz';
        end
      case 'M'
        filter='180Hz';
      case 'L'
        filter='10Hz';
      otherwise
        error(['Unknown filter char for V: ' ibvarstr(vtlen)]);
    end
  otherwise
    error(['Unknown leading char in V: ' ibvarstr(vtlen)]);
end
end