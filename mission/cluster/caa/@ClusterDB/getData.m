function out_data = getData(cdb,start_time,dt,cl_id,quantity,varargin)
%GETDATA(cdb) get Cluster data from the database or disk
% data = getData(cdb,start_epoch,dt,cl_id,quantity,options)
%
% Input:
%	cdb - ClusterDB object
%	start_epoch - start time (isdat epoch)
%	dt - time interval in sec
%	cl_id - SC#
%	quantity - one of the following:
%
%	//// EFW ////
%	e    : wE{cl_id}p12,34 -> mER
%			// electric fields (HX)
%	p    : P{cl_id} -> mPR
%			// probe potential (LX)
%	tmode: mTMode{cl_id} -> mEFWR
%			// EFW tape mode
%	dsc	 : DSC{cl_id} -> mEFWR
%			// EFW DSC
%	fdm  : FDM{cl_id} -> mEFWR
%			// EFW FDM
%	efwt : EFWT{cl_id} -> mEFWR
%			// EFW internal clock from DSC
%	ibias: IBIAS{cl_id}p{1..4} -> mEFWR
%			// EFW probe bias current
%
%	//// EFW internal burst////
%	eburst: wbE{cl_id}p12,34 -> mEFWburst
%			// electric fields 8kHz
%	pburst: P{180Hz,4kHz,32kHz}{cl_id}p{1..4} -> mEFWburstR
%			// probe potentials (180Hz,4kHz,32kHz)
%   bscburst : wBSC4kHz{cl_id} -> mEFWburstR
%			// EFW probe bias current
%
%	//// Non-standard operations ////
%   nsops:  NSOPS{cl_id} -> mEFW
%
%	//// Ephemeris ////
%	sax : SAX{cl_id} ->mEPH
%			// spin axis vector [GSE]
%	a   : A{cl_id} -> mA	// SC phase
%	r   : R{cl_id} -> mR	// SC position
%	v   : V{cl_id}, diV{cl_id} -> mR	// SC velocity
%   caa_int: write .caa_sh_interval/.caa_ms_interval file.
%
%	//// Other instruments ////
%	b   : BPP{cl_id},diBPP{cl_id}	->mBPP	// B FGM PP [GSE+DSI]
%	bfgm: B{cl_id},diB{cl_id}	->mB	// B FGM** [GSE+DSI]
%		** contact Stephan Buchert
%	bsc : wBSC{cl_id} -> mBSCR // STAFF SC B
%	edi : iEDI{cl_id},idiEDI{cl_id}	->mEDI	// E EDI PP (inert frame) [GSE+DSI]
%	ncis: NC(p,h){cl_id}			->mCIS	// N CIS PP
%	tcis: T(par,perp)C(p,h){cl_id}	->mCIS	// T CIS PP
%	vcis: VC(p,h){cl_id},diVC(p,h){cl_id}  ->mCIS	// V CIS PP [GSE+DSI]
%	wbdwf: wfWBD{cl_id} -> mWBD	// WBD waveforms E/B
%	whinat - WHINAT{cl_id} -> mWHI // WHISPER natural spectrum
%
%	options - one of the following:
%	nosave : do no save on disk
%   check_caa_sh_interval : check for .caa_sh_interval/.caa_ms_interval.
%           vcis     only fetched if .caa_sh_interval found.
%           bfgm|edi only fetched if .caa_ms_interval found.
%
% Example:
%	data = getData(...
%	ClusterDB('disco:10','/data/cluster','/tmp/my_event'),...
%	toepoch([2001 02 13 18 20 00]),120,3,'b');
%
%	This will fetch 120 sec of B FGM PP for Cluster 3 starting from
%	2001-02-13 18:20:00, using ISDAT database disco:10 or
%	CDF files in /data/cluster.
%
% See also C_GET, TOEPOCH

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(5,15)

warning off  'ISDAT:serverWarning'
warning off  'ISDAT:serverMessage'

out_data = '';

% default options
flag_save = 1;
check_caa_sh_interval=0;

for i=1:length(varargin)
  switch(varargin{i})
    case 'nosave'
      flag_save = 0;
    case 'check_caa_sh_interval'
      check_caa_sh_interval = 1;
    otherwise
      irf_log('fcal',['Option ''' varargin{i} '''not recognized'])
  end
end

save_list = '';

old_pwd = pwd;

if flag_save
  %Create the storage directory if it does not exist
  if ~exist(cdb.sp, 'dir')
    [SUCCESS,MESSAGE] = mkdir(cdb.sp);
    if SUCCESS, irf_log('save',['Created storage directory ' cdb.sp])
    else, error(MESSAGE)
    end
  end
  
  cd(cdb.sp) %enter the storage directory
  irf_log('save',[quantity ' Storage directory is ' cdb.sp])
  
  % Create .interval
  if ~exist('./.interval','file')
    fid = fopen('.interval','w');
    if fid<0, irf_log('save','problem creating .interval'),cd(old_pwd),return, end
    count = fprintf(fid,'%s %s',epoch2iso(start_time),num2str(dt));	fclose(fid);
    if count<=0, irf_log('save','problem writing to .interval'),cd(old_pwd),return, end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caa_int - write .caa_sh_interval/.caa_ms_interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'caa_int')
  load mPlan
  [iso_t,dt] = caa_read_interval;
  v_s = ['MPauseY' iso_t(1:4)];
  if ~exist(v_s,'var'), error(['Cannot load ' v_s]), end
  eval([ 'MP=' v_s ';'])
  st = iso2epoch(iso_t);
  et = st +dt;
  if any(  MP(:,1)>=st & MP(:,1)<et ) || ...
      any(  MP(:,2)>st & MP(:,2)<=et ) || ...
      any(  MP(:,1)<=st & MP(:,2)>=et ) %#ok<NODEF>
    % Create .caa_sh_interval
    fid = fopen('.caa_sh_interval','w');
    if fid<0, error('**** Problem creating .caa_sh_interval'), end
    count = fprintf(fid,'%s',epoch2iso(date2epoch(now)));
    fclose(fid);
    if count<=0,error('**** Problem writing to .caa_sh_interval'), end
    irf_log('proc','Wrote .caa_sh_interval.')
  else
    % Create .caa_ms_interval
    fid = fopen('.caa_ms_interval','w');
    if fid<0, error('**** Problem creating .caa_ms_interval'), end
    count = fprintf(fid,'%s',epoch2iso(date2epoch(now)));
    fclose(fid);
    if count<=0,error('**** Problem writing to .caa_ms_interval'), end
    irf_log('proc','Wrote .caa_ms_interval.')
  end
  flag_save=0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % nsops - check nonstandard operations table
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'nsops')
  save_file = './mEFW.mat';
  
  % Read list of nonstandard operations and see if we have one of those
  % during the requested period. Permanent problems (as loss of
  % probes, filters, etc.) must be programmed separately
  ns_ops = c_ctl('get',cl_id,'ns_ops');
  if isempty(ns_ops)
    c_ctl('load_ns_ops',[cdb.dp '/caa-control'])
    ns_ops = c_ctl('get',cl_id,'ns_ops');
  end
  if isempty(ns_ops)
    irf_log('dsrc','Nonstandard operations table not found!')
    out_data = []; cd(old_pwd), return
  end
  
  bad_start=[];
  bad_stop=[];
  bad_opcode=[];
  
  % Problem covers the whole interval
  ii = find( ns_ops(:,1)<=start_time & ns_ops(:,1)+ns_ops(:,2)>=start_time+dt );
  if ~isempty(ii)
    for j=1:length(ii)
      % no/bad data - remove the interval
      irf_log('proc',prob_s(ns_ops(ii(j),:)))
      irf_log('proc',	'whole interval in covered by ns_ops')
    end
    bad_start  =[bad_start  start_time*ones(1,length(ii))];
    bad_stop   =[bad_stop   (start_time+dt)*ones(1,length(ii))];
    bad_opcode =[bad_opcode ns_ops(ii(j),4)*ones(1,length(ii))];
  end
  % clear already processed records
  ns_ops(ii,:) = [];
  
  % Problem starts inside the interval and ends after the interval
  while 1
    ii = find( ns_ops(:,1)<start_time+dt & ns_ops(:,1)>start_time & ns_ops(:,1)+ns_ops(:,2)>=start_time+dt);
    if isempty(ii), break, end
    irf_log('proc',prob_s(ns_ops(ii(1),:)))
    irf_log('proc',	'ns_ops problem starts inside the interval and ends after the interval' )
    bad_start=[bad_start ns_ops(ii(1),1)]; %#ok<AGROW>
    bad_stop =[bad_stop  start_time+dt]; %#ok<AGROW>
    bad_opcode =[bad_opcode ns_ops(ii(1),4)]; %#ok<AGROW>
    % clear already processed records
    ns_ops(ii(1),:) = [];
  end
  
  % Problem starts before the interval and ends inside the interval
  while 1
    ii = find( ns_ops(:,1)<=start_time & ns_ops(:,1)+ns_ops(:,2)>start_time & ns_ops(:,1)+ns_ops(:,2)<=start_time+dt);
    if isempty(ii), break, end
    irf_log('proc',prob_s(ns_ops(ii(1),:)))
    irf_log('proc',	'ns_ops problem starts before the interval and ends inside the interval' )
    bad_start=[bad_start start_time]; %#ok<AGROW>
    bad_stop =[bad_stop  ns_ops(ii(1),1)+ns_ops(ii(1),2)]; %#ok<AGROW>
    bad_opcode =[bad_opcode ns_ops(ii(1),4)]; %#ok<AGROW>
    ns_ops(ii(1),:) = [];
  end
  
  % Problem is inside the interval
  while 1
    ii = find( ns_ops(:,1)>start_time & ns_ops(:,1)+ns_ops(:,2)<start_time+dt);
    if isempty(ii), break, end
    irf_log('proc',prob_s(ns_ops(ii(1),:)))
    irf_log('proc',	'ns_ops problem inside the interval' )
    bad_start=[bad_start ns_ops(ii(1),1)]; %#ok<AGROW>
    bad_stop =[bad_stop  ns_ops(ii(1),1)+ns_ops(ii(1),2)]; %#ok<AGROW>
    bad_opcode =[bad_opcode ns_ops(ii(1),4)]; %#ok<AGROW>
    ns_ops(ii(1),:) = [];
  end
  
  bad_intervals = [bad_start' bad_stop' bad_opcode']; %#ok<NASGU>
  c_eval('NSOPS?=bad_intervals;save_list=[save_list '' NSOPS? ''];',cl_id);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dsc - EFW DSC old
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'dscold')
  save_file = './mEFWR.mat';
  
  [t,dsc] = caa_is_get(cdb.db, start_time, dt, cl_id, 'efw', 'DSC');
  if isempty(dsc)
    irf_log('dsrc',irf_ssub('No data for DSC?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  
  % DSC fields we want to save
  % 0:55 - format tables
  % 64:68 - HX, LX format pointers
  % 73 - executive version
  % 78 - sc #
  % 79 - trap_cointer
  % 80:84 - clock
  % 128:139 - bias, stub, guard settings
  % 210:249 - sweep settings
  dsc_is = [0:55 64:68 73 78 79 80:84 128:139 210:249] + 1;
  dsc_i = [0:55 64:68 73 78 79 128:139] + 1;
  
  % storage variables
  t_start_save = [];
  dsc_save = [];
  jump_flag = [];
  n_good = 0;
  n_jumpy = 0;
  
  % temporal variables
  t_st = [];
  t_end = [];
  dsc_good = [];
  
  dsc_last = dsc(:,1);
  t_dsc_last = t(1);
  count_good = -1;
  
  for i = 1:length(t)
    if sum(abs( dsc(dsc_i,i) - dsc_last(dsc_i))) == 0
      % Same as previous
      count_good = count_good + 1;
      if count_good == 1
        t_st = t_dsc_last;
        dsc_good = dsc_last;
      end
      t_end = t(i);
    else
      % Differs from previous
      dsc_last = dsc(:,i);
      l = length(jump_flag);
      if count_good == 0
        % The previous point was also different
        t_start_save(l+1) = t(i); %#ok<AGROW>
        dsc_save(:,l+1) = dsc(dsc_is,i); %#ok<AGROW>
        jump_flag(l+1) = 1; %#ok<AGROW>
        n_jumpy = n_jumpy + 1;
      else
        % Save a good interval
        t_start_save(l+1) = t_st; %#ok<AGROW>
        dsc_save(:,l+1) = dsc_good(dsc_is); %#ok<AGROW>
        jump_flag(l+1) = 0; %#ok<AGROW>
        n_good = n_good + 1;
        irf_log('dsrc',['Saving good from ' ...
          epoch2iso(t_st,1) '-' epoch2iso(t_end,1)])
      end
      count_good = 0;
    end
    t_dsc_last = t(i);
  end
  
  if count_good > 0 % The last interval was also good
    % Save a good interval
    t_start_save(end+1) = t_st; %#ok<NASGU>
    dsc_save(:,end+1) = dsc_good(dsc_is); %#ok<NASGU>
    n_good = n_good + 1;
    irf_log('dsrc',['Saving good from ' epoch2iso(t_st,1) '-' epoch2iso(t_end,1)])
  end
  
  irf_log('dsrc',sprintf('\nFound total %d good and %d jumpy intervals',...
    n_good, n_jumpy))
  c_eval('DSC?=[t_start_save'' dsc_save'']'';save_list=[save_list '' DSC? ''];',cl_id);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dsc - EFW DSC new save all
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'dsc')
  save_file = './mEFWR.mat';
  
  [t,dsc] = caa_is_get(cdb.db, start_time, dt, cl_id, 'efw', 'DSC');
  if isempty(t) || isempty(dsc)
    irf_log('dsrc',irf_ssub('No data for DSC?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  
  % All DSC fields are saved
  
  c_eval('DSC?=[t dsc''];save_list=[save_list '' DSC? ''];',cl_id);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % fdm - EFW FDM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'fdm')
  save_file = './mEFWR.mat';
  
  [t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, 'efw', 'FDM'); %#ok<ASGLU>
  if isempty(data)
    irf_log('dsrc',irf_ssub('No data for FDM?',cl_id))
    out_data = []; cd(old_pwd), return
  else, c_eval('FDM?=[t data''];',cl_id);
  end
  
  c_eval('save_list=[save_list '' FDM? ''];',cl_id);
  clear t data
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ibias - EFW probe bias current
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'ibias')
  save_file = './mEFWR.mat';
  
  probe_list = 1:4;
  p_ok = [];
  
  % Check for p1 problems on SC2,3
  if (start_time>toepoch([2002 07 29 09 06 59 ]) && cl_id==3) || ...
      (start_time>toepoch([2007 05 13 03 23 30]) && cl_id==2)
    probe_list = 2:4;
    irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
  end
  
  % Check for p1, p4 problems on SC1
  if cl_id==1
    if start_time > toepoch([2009 10 14 07 00 00]) || ...
        (start_time > toepoch([2009 04 19 00 00 00]) && start_time < toepoch([2009 05 07 00 00 00]))
      probe_list = 2:3;
      irf_log('dsrc',sprintf('p1 and p4 are BAD on sc%d',cl_id));
    elseif start_time>toepoch([2001 12 28 03 00 00])
      probe_list = 2:4;
      irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
    end
  end
  
  for probe=probe_list
    [t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, ...
      'efw', 'E', ['p' num2str(probe)],'bias'); %#ok<ASGLU>
    if isempty(data)
      irf_log('dsrc',irf_ssub('No data for IBIAS?p!',cl_id,probe))
    else
      eval(irf_ssub('IBIAS?p!=[t data];',cl_id,probe))
      p_ok = [p_ok probe]; %#ok<AGROW>
    end
    clear t data
  end
  
  if isempty(p_ok), out_data = []; cd(old_pwd), return, end
  for probe=p_ok
    eval(irf_ssub('save_list=[save_list ''IBIAS?p! ''];',cl_id,probe))
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % efwt - EFW clock
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'efwt')
  save_file = './mEFWR.mat';
  
  % Read EFW clock to check for time since last reset
  t = []; out_data = []; efwtime = [];
  for st_tmp = start_time-16:32:start_time+dt+16
    [t_tmp,data] = caa_is_get(cdb.db, st_tmp, 32, cl_id, 'efw', 'DSC');
    if ~isempty(data)
      efwtime_tmp = (data(81,:) +data(82,:)*256 +data(83,:)*65536 + ...
        data(84,:)*16777216 +data(85,:)*4294967296)/1000;
      efwtime = [efwtime efwtime_tmp]; %#ok<AGROW>
      t = [t t_tmp']; %#ok<AGROW>
    end
  end
  
  if isempty(efwtime)
    irf_log('dsrc',irf_ssub('No data for EFWT?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  c_eval(['EFWT?=[t; efwtime]'';'...
    'save_list=[save_list '' EFWT? ''];'],cl_id);
  clear t data efwtime
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % tmode - EFW tape mode
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'tmode')
  save_file = './mEFWR.mat';
  
  % Find TapeMode
  % We read FDM from isdat and 5-th column contains the HX mode
  % (undocumented feature)
  % 0 - normal mode  (V12L,V34L)
  % 1 - tape mode 1  (V12M,V34M)
  % 2 - tape mode 2  (V12M,V34M)
  % 3 - tape mode 3  (V1M,V2M,V3M,V4M)
  [t,data] = caa_is_get(cdb.db,start_time,dt,cl_id,'efw','FDM'); %#ok<ASGLU>
  if isempty(data)
    irf_log('dsrc',irf_ssub('No data for mTMode?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  
  c_eval(['mTMode?=data(5,:);'...
    'save_list=[save_list '' mTMode? ''];'],cl_id)
  clear t data
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % e - Electric field
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'e') || strcmp(quantity,'eburst')
  
  if strcmp(quantity,'eburst'), do_burst = 1; else, do_burst = 0; end
  if do_burst
    save_file = './mEFWburstR.mat';
    tmmode='burst';
    param='8kHz';
    var_name = 'wbE?p';
  else
    save_file = './mER.mat';
    tmmode='hx';
    var_name = 'wE?p';
    
    [ok,tm] = c_load('mTMode?',cl_id);
    if ~ok
      tmode = getData(cdb,start_time,dt,cl_id,'tmode');
      if isempty(tmode)
        irf_log('dsrc',irf_ssub('Cannot load mTMode?',cl_id))
        out_data = []; cd(old_pwd), return
      end
      tm = tmode{2};
      clear tmode
    end
    if isempty(tm)
      irf_log('dsrc',irf_ssub('Cannot load mTMode?',cl_id))
      out_data = []; cd(old_pwd), return
    end
    
    if any(tm~=tm(1))
      irf_log('dsrc','tape mode changes during the selected time inteval')
      irf_log('dsrc','data interval will be truncated')
    end
    tm = tm(1);
    if tm<1e-30, param='10Hz'; else, param='180Hz'; end
    clear tm
    
    %%%%%%%%%%%%%%%%%%%%%%%%% FILTER MAGIC %%%%%%%%%%%%%%%%%%%%%%
    switch cl_id
      case 1
        if start_time>toepoch([2015 02 26 09 35 00])
          param='180Hz';
        end
      case 2
        if start_time>toepoch([2001 07 23 13 54 18])
          % 10Hz filter problem on SC2
          param='180Hz';
        elseif cl_id==2 && start_time<toepoch([2001 07 23 13 54 18]) && ...
            start_time+dt>toepoch([2001 07 23 13 54 18])
          % Request interval overlaps with the time when the 10Hz filter
          % got broken on SC2. We truncate the request.
          dt = toepoch([2001 07 23 13 54 17]) - start_time;
          for in=1:length(start_time)
            if start_time(in)<toepoch([2001 07 23 13 54 18]) && ...
                start_time(in)+dt(in)>toepoch([2001 07 23 13 54 18])
              dt(in) = toepoch([2001 07 23 13 54 17]) ...
                - start_time(in);
            end
          end
          irf_log('proc', ...
            '10Hz filter got broken inside the requested interval')
          irf_log('proc',	['truncating interval: setting DT to ' num2str(dt)])
        elseif (((cl_id==1 && start_time>toepoch([2001 07 30 17 05 54.9])) || ...
            (cl_id==3 && start_time>toepoch([2001 07 31 00 12 29.5]))) && ...
            start_time<toepoch([2001 09 02 23 15 00])) || ...
            (cl_id==4 && ((start_time>toepoch([2001 07 31 04 55 33.15]) && ...
            start_time<toepoch([2001 08 02 11 25 40])) || ...
            (start_time>toepoch([2001 08 06 23 58 50.7]) && ...
            start_time<toepoch([2001 09 02 23 15 00]))))
          % all sc run on 180Hz filter in august 2001 most of the time
          param='180Hz';
        elseif start_time>toepoch([2001 09 10 04 21 57.6]) && ...
            start_time<toepoch([2001 09 17 05 27 54])
          % this needs to be investigated....
          param='180Hz';
        end
      case 4
        if start_time>toepoch([2015 02 28 13 00 00])
          param='180Hz';
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%% END FILTER MAGIC %%%%%%%%%%%%%%%%%%%%
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%% PROBE MAGIC %%%%%%%%%%%%%%%%%%%%%%
  pl = [12,34];
  switch cl_id
    case 1
      if start_time>toepoch([2018 12 10 03 00 16])
        pl = [];
        irf_log('dsrc',sprintf('            !No diff measurement on sc%d',cl_id));
      elseif start_time>toepoch([2009 10 14 07 00 00]) ||  ...
          (start_time>toepoch([2009 04 19 00 00 00]) && start_time<toepoch([2009 05 07 00 00 00]))
        pl = 32;
        irf_log('dsrc',sprintf('  !Only p32 exists on sc%d',cl_id));
      elseif (start_time>toepoch([2003 9 29 00 27 0]) || ...
          (start_time>toepoch([2003 3 27 03 50 0]) && start_time<toepoch([2003 3 28 04 55 0])) ||...
          (start_time>toepoch([2003 4 08 01 25 0]) && start_time<toepoch([2003 4 09 02 25 0])) ||...
          (start_time>toepoch([2003 5 25 15 25 0]) && start_time<toepoch([2003 6 08 22 10 0])) )
        pl = [32, 34];
        irf_log('dsrc',sprintf('  !Using p32 on sc%d',cl_id));
      elseif start_time>toepoch([2001 12 28 03 00 00])
        pl = 34;
        irf_log('dsrc',sprintf('  !Only p34 exists on sc%d',cl_id));
      elseif  (start_time>=toepoch([2001 04 12 03 00 00]) && start_time<toepoch([2001 04 12 06 00 00])) || ...
          (  start_time>=toepoch([2001 04 14 06 00 00]) && start_time<toepoch([2001 04 16 15 00 00])) || ...
          (  start_time>=toepoch([2001 04 18 03 00 00]) && start_time<toepoch([2001 04 20 09 00 00])) || ...
          (  start_time>=toepoch([2001 04 21 21 00 00]) && start_time<toepoch([2001 04 22 03 00 00])) || ...
          (  start_time>=toepoch([2001 04 23 09 00 00]) && start_time<toepoch([2001 04 23 15 00 00]))
        % The bias current is a bit too large
        % on p3 and p4 on C1&2 in April 2001.
        % Ignore p3, p4 and p34 and only use p1, p2 and p12.
        % Use only complete 3-hour intervals to keep it simple.
        pl = 12;
        irf_log('dsrc',sprintf('            !Too high bias current on p34 for sc%d',cl_id));
      end
    case 2
      if start_time>toepoch([2015 10 12 12 00 0])
        pl = 34;
        irf_log('dsrc',sprintf('  !Only p34 exists on sc%d',cl_id));
      elseif start_time>toepoch([2007 11 24 15 40 0])
        pl = [32, 34];
        irf_log('dsrc',sprintf('  !Using p32 on sc%d',cl_id));
      elseif start_time>toepoch([2007 05 13 03 23 48])
        pl = 34;
        irf_log('dsrc',sprintf('  !Only p34 exists on sc%d',cl_id));
      elseif (start_time>=toepoch([2001 04 09 21 00 00]) && start_time<toepoch([2001 04 10 06 00 00])) || ...
          (  start_time>=toepoch([2001 04 10 09 00 00]) && start_time<toepoch([2001 04 19 15 00 00])) || ...
          (  start_time>=toepoch([2001 04 20 03 00 00]) && start_time<toepoch([2001 04 23 15 00 00])) || ...
          (  start_time>=toepoch([2001 04 24 00 00 00]) && start_time<toepoch([2001 04 24 15 00 00]))
        pl = 12;
        irf_log('dsrc',sprintf('  !Too high bias current on p34 for sc%d',cl_id));
      end
    case 3
      if start_time>toepoch([2011 6 01 09 30 0])
        pl = [];
        irf_log('dsrc',sprintf('            !No diff measurement on sc%d',cl_id));
      elseif start_time>toepoch([2003 9 29 00 27 0]) || ...
          (start_time>toepoch([2003 3 27 03 50 0]) && start_time<toepoch([2003 3 28 04 55 0])) ||...
          (start_time>toepoch([2003 4 08 01 25 0]) && start_time<toepoch([2003 4 09 02 25 0])) ||...
          (start_time>toepoch([2003 5 25 15 25 0]) && start_time<toepoch([2003 6 08 22 10 0]))
        pl = [32, 34];
        irf_log('dsrc',sprintf('  !Using p32 on sc%d',cl_id));
      elseif start_time>toepoch([2002 07 29 09 06 59])
        pl = 34;
        irf_log('dsrc',sprintf('  !Only p34 exists on sc%d',cl_id));
      end
    case 4
      if start_time>=toepoch([2013 07 01 13 30 00]) % 2013 07 01 14 43 44
        pl = 12;
        irf_log('dsrc',sprintf('  !Only p12 exists on sc%d',cl_id));
      end
  end
  if isempty(pl), out_data = []; cd(old_pwd), return, end
  
  %%%%%%%%%%%%%%%%%%%%%%% END PROBE MAGIC %%%%%%%%%%%%%%%%%%%%
  
  p_ok = [];
  for probe=pl
    irf_log('dsrc',['EFW...sc' num2str(cl_id)...
      '...Ep' num2str(probe) ' ' param ' filter']);
    t = [];	data = [];
    for in=1:length(start_time)
      if length(start_time)>1
        irf_log('dsrc',...
          sprintf('chunk #%d : %s %d sec',in,...
          epoch2iso(start_time(in),1),dt(in)))
      end
      [t_tmp,data_tmp] = caa_is_get(cdb.db, start_time(in), dt(in), cl_id, ...
        'efw', 'E', ['p' num2str(probe)], param, tmmode);
      t = [t; t_tmp]; data = [data; data_tmp]; clear t_tmp data_tmp %#ok<AGROW>
    end
    if ~isempty(data)
      % Correct start time of the burst
      if do_burst
        burst_f_name = [cdb.dp '/burst/'...
          irf_ssub([irf_fname(t(1),1) 'we.0?'],cl_id)];
        if exist(burst_f_name,'file')
          err_t = t(1) - c_efw_burst_chkt(cdb.db,burst_f_name);
          
          t = t - err_t;
          irf_log('dsrc',...
            ['burst start time was corrected by ' num2str(err_t) ' sec'])
        else, irf_log('dsrc','burst start time was not corrected')
        end
      end
      
      data = check_timeline([t data]); %#ok<NASGU>
      c_eval([var_name  num2str(probe) '=data;'],cl_id);
      p_ok = [p_ok probe]; %#ok<AGROW>
      
    else
      irf_log('dsrc',...
        irf_ssub(['No data for ' var_name num2str(probe)],cl_id))
    end
  end
  
  if isempty(p_ok), out_data = []; cd(old_pwd), return, end
  for probe=p_ok
    c_eval(['save_list=[save_list '' ' var_name num2str(probe) '''];'],cl_id);
  end
  clear t data tm pl param
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % p - SC potential
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'p') || strcmp(quantity,'pburst')
  
  if strcmp(quantity,'pburst'), do_burst = 1; else, do_burst = 0; end
  if do_burst
    save_file = './mEFWburstR.mat';
    tmmode='burst';
    param={'180Hz','4kHz','32kHz'};
  else
    save_file = './mPR.mat';
    param={'10Hz'}; tmmode='lx';
  end
  
  probe_list = 1:4;
  
  %%%%%%%%%%%%%%%%%%%%%%%%% PROBE MAGIC %%%%%%%%%%%%%%%%%%%%%%
  switch cl_id
    case 1
      if start_time>toepoch([2018 12 10 03 00 16])
        % p2 failure
        probe_list = 2;
        irf_log('dsrc',sprintf('p1, p3 and p4 are BAD on sc%d',cl_id))
      elseif start_time>toepoch([2009 10 14 07 00 00]) || ...
          (start_time>toepoch([2009 04 19 00 00 00]) && start_time<toepoch([2009 05 07 00 00 00]))
        % p4 failure
        probe_list = 2:3;
        irf_log('dsrc',sprintf('p1 and p4 are BAD on sc%d',cl_id))
      elseif start_time>toepoch([2001 12 28 03 00 00])
        % p1 failure
        probe_list = 2:4;
        irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
      elseif ( (start_time>=toepoch([2001 04 12 03 00 00]) && start_time<toepoch([2001 04 12 06 00 00])) || ...
          (  start_time>=toepoch([2001 04 14 06 00 00]) && start_time<toepoch([2001 04 16 15 00 00])) || ...
          (  start_time>=toepoch([2001 04 18 03 00 00]) && start_time<toepoch([2001 04 20 09 00 00])) || ...
          (  start_time>=toepoch([2001 04 21 21 00 00]) && start_time<toepoch([2001 04 22 03 00 00])) || ...
          (  start_time>=toepoch([2001 04 23 09 00 00]) && start_time<toepoch([2001 04 23 15 00 00])) )
        % The bias current is a bit too large
        % on p3 and p4 on C1&2 in April 2001.
        % Ignore p3, p4 and p34 and only use p1, p2 and p12.
        % Use only complete 3-hour intervals to keep it simple.
        probe_list = [1 2];
        irf_log('dsrc',sprintf('Too high bias current on p3&p4 sc%d',cl_id));
      end
      if start_time>toepoch([2015 02 26 09 35 00])
        param={'180Hz'};
      end
    case 2
      if start_time>toepoch([2015 10 12 12 00 0])
        % P2 failure
        probe_list = 3:4;
        irf_log('dsrc',sprintf('p1&p2 are BAD on sc%d',cl_id))
      elseif start_time>=toepoch([2007 06 01 17 20 00])
        % We use 180 Hz filter
        if ~do_burst, param={'180Hz'}; end
        irf_log('dsrc',sprintf('using 180Hz filter on sc%d',cl_id))
        probe_list = 2:4;
        irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
      elseif start_time>=toepoch([2007 05 13 03 23 48])
        probe_list = [2 4];
        irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id))
        irf_log('dsrc',sprintf('10Hz filter problem on p3 sc%d',cl_id))
      elseif start_time+dt>toepoch([2001 07 23 13 54 18]) && ~do_burst
        % 10Hz filter problem on C2 p3
        % Any changes should also go to ClusterProc/getData/probesa
        probe_list = [1 2 4];
        irf_log('dsrc',sprintf('10Hz filter problem on p3 sc%d',cl_id))
      elseif ( (start_time>=toepoch([2001 04 09 21 00 00]) && start_time<toepoch([2001 04 10 06 00 00])) || ...
          (  start_time>=toepoch([2001 04 10 09 00 00]) && start_time<toepoch([2001 04 19 15 00 00])) || ...
          (  start_time>=toepoch([2001 04 20 03 00 00]) && start_time<toepoch([2001 04 23 15 00 00])) || ...
          (  start_time>=toepoch([2001 04 24 00 00 00]) && start_time<toepoch([2001 04 24 15 00 00])) )
        % The bias current is a bit too large
        % on p3 and p4 on C1&2 in April 2001.
        % Ignore p3, p4 and p34 and only use p1, p2 and p12.
        % Use only complete 3-hour intervals to keep it simple.
        probe_list = [1 2];
        irf_log('dsrc',sprintf('Too high bias current on p3&p4 sc%d',cl_id));
      end
      if start_time>toepoch([2015 02 26 09 35 00])
        param={'180Hz'};
      end
    case 3
      if start_time>toepoch([2014 11 03 20 58 16.7])
        % p2 failure
        probe_list = 4;
        irf_log('dsrc',sprintf('p1, p2 & p3 are BAD on sc%d',cl_id));
      elseif start_time>toepoch([2011 6 01 09 30 0])
        % p3 failure
        probe_list = [2 4];
        irf_log('dsrc',sprintf('p1 & p3 are BAD on sc%d',cl_id));
      elseif start_time>toepoch([2002 07 29 09 06 59 ])
        % p1 failure
        probe_list = 2:4;
        irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
      end
      if start_time>toepoch([2015 03 08 04 10 00])
        param={'180Hz'};
      end
    case 4
      if start_time>=toepoch([2015 02 17 07 30 00]) % 2015-02-17 07:36:30
        probe_list = 1:2;
        irf_log('dsrc',sprintf('p3 & p4 are BAD on sc%d',cl_id));
      elseif start_time>=toepoch([2013 07 01 13 30 00]) % 2013 07 01 14 43 44
        probe_list = 1:3;
        irf_log('dsrc',sprintf('p4 is BAD on sc%d',cl_id));
      end
      if start_time>toepoch([2015 03 08 04 10 00])
        param={'180Hz'};
      end
  end
  %%%%%%%%%%%%%%%%%%%%%%% END PROBE MAGIC %%%%%%%%%%%%%%%%%%%%
  
  n_ok = 0;
  for j=1:length(param)
    for probe=probe_list
      irf_log('dsrc',['EFW...sc' num2str(cl_id) '...probe' num2str(probe)...
        '->P' param{j} num2str(cl_id) 'p' num2str(probe)]);
      t = [];	data = [];
      for in=1:length(start_time)
        if length(start_time)>1
          irf_log('dsrc',...
            sprintf('chunk #%d : %s %d sec',in,...
            epoch2iso(start_time(in),1),dt(in)))
        end
        [t_tmp,data_tmp] = caa_is_get(cdb.db, start_time(in), dt(in), cl_id, ...
          'efw', 'E', ['p' num2str(probe)],param{j}, tmmode);
        t = [t; t_tmp]; data = [data; data_tmp]; clear t_tmp data_tmp %#ok<AGROW>
      end
      
      if ~isempty(data)
        % Correct start time of the burst
        if do_burst
          burst_f_name = [cdb.dp '/burst/'...
            irf_ssub([irf_fname(t(1),1) 'we.0?'],cl_id)];
          if exist(burst_f_name,'file')
            start_satt = c_efw_burst_chkt(cdb.db,burst_f_name);
            if isempty(start_satt)
              irf_log('dsrc','burst start time was not corrected')
            else
              err_t = t(1) - start_satt;
              irf_log('dsrc',['burst start time was corrected by ' ...
                num2str(err_t) ' sec'])
              t = t - err_t;
            end
          else
            irf_log('dsrc','burst start time was not corrected')
          end
        end
        
        data = check_timeline([t data]); %#ok<NASGU>
        if ~do_burst
          % XXX: Hack to avoid saving variables named P180HzXpXX
          % which are sampled on SC2 since June 1, 2007
          % A lot of software depends on LX data to be called
          % P10HzXpXX
          eval(irf_ssub(['P10Hz?p!=data;'...
            'save_list=[save_list ''P10Hz?p! ''];'],cl_id,probe));
        else
          eval(irf_ssub(['P' param{j} '?p!=data;'...
            'save_list=[save_list ''P' param{j} '?p! ''];'],cl_id,probe));
        end
        n_ok = n_ok + 1;
        
      else, irf_log('dsrc', irf_ssub(['No data for P' param{j} '?p!'],cl_id,probe));
      end
      clear t data
    end
  end
  
  if ~n_ok, out_data = []; cd(old_pwd), return, end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % aux data - Phase, etc.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'a')
  save_file = './mA.mat';
  pha = [];
  downloadStatus = 0;
  irf_log('dsrc','Trying to to read phase from CP_AUX_SPIN_TIME...');
  currentDir = pwd;	tempDir = sprintf('CAA_Download_%d',fix(rand*1e6));
  mkdir(tempDir); cd(tempDir);
  tint = start_time +[-5 dt+10];
  datasetName = sprintf('C%d_CP_AUX_SPIN_TIME',cl_id);
  try
    downloadStatus = caa_download(tint,datasetName,'stream');
  catch, irf_log('dsrc','Error streaming from CSA')
  end
  d = dir(['CAA/' datasetName '/*.cef.gz']);
  if ~isempty(d) && downloadStatus
    cefFile = ['CAA/' datasetName '/' d.name];
    cef_init(); cef_read(cefFile);
    c1 = onCleanup(@() cef_close());
    c2 = onCleanup(@() rmdir([currentDir '/' tempDir],'s'));
    tt = cef_var('time_tags');
    if isempty(tt) % check for empty cef file return
      irf_log('dsrc','did not suceed: zero data points returned')
    else
      tt = irf_time( cef_date(tt'),'datenum>epoch');
      spinPeriod = cef_var('spin_period'); spinPeriod = double(spinPeriod');
      % find errors
      iJump = find(abs(spinPeriod-median(spinPeriod))>5*std(spinPeriod));
      if length(iJump) < min(4,length( spinPeriod ))
        if length(iJump)>1 && iJump(end)-iJump(1)==2, iJump=iJump(1)+(1:3)'-1; end
        for i = iJump'
          irf_log('proc',['removing erroneous point at ' epoch2iso(tt(i))])
        end
        spinPeriod(iJump) = [];
        tt(iJump) = [];
      end
      refTime = [-.5 .5]; % part of spin
      refPhase = refTime*360+180; % spin period center corresponds to phase 0
      deltaT = repmat(refTime,size(tt,1),1).*...
        repmat(double(spinPeriod),1,length(refTime));
      tmat = repmat(tt(:,1),1,length(refTime))+deltaT;
      amat = repmat(refPhase,size(tt,1),1);
      tmat = reshape(tmat',numel(tmat),1);
      difftmat = diff(tmat); ii = find(difftmat<0); tmat(ii) = tmat(ii+1);
      if sum(tmat>=start_time & tmat<=start_time+dt)>1 % at least 2 points
        amat = reshape(amat',numel(amat),1);
        pha = [tmat amat];
      else
        irf_log('dsrc','did not suceed: too few data points returned')
      end
    end
  end
  cd(currentDir)
  if isempty(pha) % read from isdat
    irf_log('dsrc','Reading phase from ISDAT instead');
    % We ask for 2 sec more from each side
    % to avoid problems with interpolation.
    [t,data] = c_get_phase(cdb.db,start_time-2,dt+4,cl_id,'phase_2');
    if ~isempty(data) && length(t)>1, pha=[t data]; end
    clear t data
  end
  if ~isempty(pha)
    c_eval('Atwo?=pha;save_list=[save_list '' Atwo? ''];',cl_id);
  else
    out_data = [];
    c_eval('Atwo?=[];save_list=[save_list '' Atwo? ''];',cl_id);
    irf_log('dsrc',irf_ssub('No/short data for Atwo?',cl_id))
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % aux data - Position
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'r')
  save_file = './mR.mat';
  
  [t,data] = caa_is_get(cdb.db, start_time, dt, ...
    cl_id, 'ephemeris', 'position');
  if ~isempty(data)
    % Remove points exactly equal to the end time
    % to avoid duplicate time at the start of the next interval
    if (t(end) >= start_time+dt-0.001) && (length(t) > 1)
      t=t(1:end-1); %#ok<NASGU>
      data=data(:,1:end-1); %#ok<NASGU>
    end
    c_eval('R?=[t data''];save_list=[save_list '' R? ''];',cl_id);
  else
    irf_log('dsrc',irf_ssub('No data for R?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  clear t data
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % aux data - Velocity
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'v')
  save_file = './mR.mat';
  
  [t,data] = caa_is_get(cdb.db, start_time, dt, ...
    cl_id, 'ephemeris', 'velocity'); %#ok<ASGLU>
  if isempty(data)
    irf_log('dsrc',irf_ssub('No data for V?, diV?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  
  c_eval('V?=[t data''];save_list=[save_list '' V? ''];',cl_id);
  clear t data
  
  % Transform vector data to DSI
  [ok,sax] = c_load('SAX?',cl_id);
  if ~ok
    tempv = getData(cdb,start_time,dt,cl_id,'sax');
    if isempty(tempv)
      irf_log('dsrc',irf_ssub('Cannot load SAX?',cl_id))
      sax = [];
    else, sax = tempv{2};
    end
    clear tempv
  end
  if ~isempty(sax)
    c_eval('diV?=c_gse2dsi(V?,sax);save_list=[save_list '' diV? ''];',cl_id);
  else, irf_log('dsrc',irf_ssub('No data for diV?',cl_id))
  end
  
  %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magc - location in magnetic coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'magc')
	save_file='./mEPH.mat';
	
	var_list = {'lt', 'mlt', 'l_shell','inv_lat'};
	var_list_s = {'LT?', 'MLT?', 'L?','ILAT?'};
	
	for j=1:length(var_list)
		[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, 'ephemeris', var_list{j});
		if ~isempty(data)
			c_eval([var_list_s{j} '=[double(t) double(data'')];'],cl_id);
			clear t data;
			c_eval(['save_list=[save_list '' ' var_list_s{j} ' ''];'],cl_id);
		else
			irf_log('dsrc',irf_ssub(['No data for ' var_list_s{j}],cl_id))
		end
	end
  %}
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B FGM - full res LOCAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'bfgmlocal')
  save_file = './mB.mat';
  
  if check_caa_sh_interval
    if ~exist('./.caa_ms_interval','file')
      irf_log('proc','Outside magnetosphere. No bfgm data fetched.')
      out_data = []; cd(old_pwd), return
    end
  end
  
  dat = c_get_bfgm(start_time + [0 dt],cl_id);
  
  if isempty(dat)
    irf_log('dsrc',irf_ssub('No data for B, diB?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  c_eval('B?=dat;save_list=[save_list '' B? ''];',cl_id);
  clear dat
  
  % Transform vector data to DSI
  [ok,sax] = c_load('SAX?',cl_id);
  if ~ok
    tempv = getData(cdb,start_time,dt,cl_id,'sax');
    if isempty(tempv)
      irf_log('dsrc',irf_ssub('Cannot load SAX?',cl_id))
      sax = [];
    else, sax = tempv{2};
    end
    clear tempv
  end
  if ~isempty(sax)
    c_eval('diB?=c_gse2dsi(B?,sax);save_list=[save_list '' diB? ''];',cl_id);
  else, irf_log('dsrc',irf_ssub('No data for diB?',cl_id))
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B FGM - full res from CAA
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'bfgm')
  save_file = './mB.mat';
  
  if check_caa_sh_interval
    if ~exist('./.caa_ms_interval','file')
      irf_log('proc','Outside magnetosphere. No bfgm data fetched.')
      out_data = []; cd(old_pwd), return
    end
  end
  
  downloadStatus = 0;
  currentDir = pwd; tempDir = tempname; dat = [];
  try
    dsetName = irf_ssub('C?_CP_FGM_FULL',cl_id);
    mkdir(tempDir);
    cd(tempDir);
    downloadStatus = caa_download(start_time + [0 dt],dsetName,'stream');
    if downloadStatus
      cd(['CAA' filesep dsetName]);
      d=dir('*.cef.gz');
      dat = c_caa_cef_var_get('B_vec_xyz_gse',d.name);
    end
  catch
    irf_log('dsrc',irf_ssub('Error downloading B from CAA',cl_id))
  end
  cd(currentDir);
  if exist(tempDir,'dir'), rmdir(tempDir,'s'); end
  
  if isempty(dat)
    irf_log('dsrc',irf_ssub('No data for B, diB?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  c_eval('B?=dat;save_list=[save_list '' B? ''];',cl_id);
  clear dat
  
  % Transform vector data to DSI
  [ok,sax] = c_load('SAX?',cl_id);
  if ~ok
    tempv = getData(cdb,start_time,dt,cl_id,'sax');
    if isempty(tempv)
      irf_log('dsrc',irf_ssub('Cannot load SAX?',cl_id))
      sax = [];
    else, sax = tempv{2};
    end
    clear tempv
  end
  if ~isempty(sax)
    c_eval('diB?=c_gse2dsi(B?,sax);save_list=[save_list '' diB? ''];',cl_id);
  else, irf_log('dsrc',irf_ssub('No data for diB?',cl_id))
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B STAFF SC
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
elseif strcmp(quantity,'bsc')
  
  save_file = './mBSCR.mat';
  
  [ok,tm] = c_load('mTMode?',cl_id);
  if ~ok
    tmode = getData(cdb,start_time,dt,cl_id,'tmode');
    if isempty(tmode)
      irf_log('dsrc',irf_ssub('Cannot load mTMode?',cl_id))
      out_data = []; cd(old_pwd), return
    end
    tm = tmode{2};
    clear tmode
  end
  if isempty(tm)
    irf_log('dsrc',irf_ssub('Cannot load mTMode?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  
  if any(tm~=tm(1))
    irf_log('dsrc','tape mode changes during the selected time inteval')
    irf_log('dsrc','data interval will be truncated')
  end
  tm = tm(1);
  if tm<1e-30, param='10Hz'; else, param='180Hz'; end
  clear tm
  
  [t,data] = caa_is_get(cdb.db, start_time, dt, ...
    cl_id, 'staff', 'B_SC','Bx_By_Bz', ['0-' param]); %#ok<ASGLU>
  if isempty(data)
    irf_log('dsrc',irf_ssub('No data for wBSC?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  
  c_eval('wBSC?=[t data''];save_list=[save_list '' wBSC? ''];',cl_id);
  clear t data
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B STAFF SC in the EFW IB
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'bscburst')
  
  save_file = './mEFWburstR.mat';
  
  co = 'xyz';
  for comp=1:3
    [t,data] = caa_is_get(cdb.db, start_time, dt, cl_id,...
      'efw','dB',co(comp),'4kHz','burst','tm');
    if isempty(t) || isempty(data)
      irf_log('dsrc',irf_ssub('No data for wBSC4kHz?',cl_id))
      out_data = []; cd(old_pwd), return
    end
    B(:,comp) = data; %#ok<AGROW>
  end
  
  B = -B/7000; % Convert to V - same as STAFF B_SC_Level_1
  
  % Correct start time of the burst
  burst_f_name = [cdb.dp '/burst/'...
    irf_ssub([irf_fname(t(1),1) 'we.0?'],cl_id)];
  if exist(burst_f_name,'file')
    err_t = t(1) - c_efw_burst_chkt(cdb.db,burst_f_name);
    
    t = t - err_t;
    irf_log('dsrc',...
      ['burst start time was corrected by ' num2str(err_t) ' sec'])
  else, irf_log('dsrc','burst start time was not corrected')
  end
  
  B = rm_ib_spike([t B]);
  
  B = c_efw_burst_bsc_tf(B,cl_id); %#ok<NASGU> % Apply the transfer function
  
  c_eval('wBSC4kHz?=B;save_list=[save_list '' wBSC4kHz? ''];',cl_id);
  clear t data
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CSDS PP [GSE+DSI]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'b') || strcmp(quantity,'edi') || ...
    strcmp(quantity,'ncis') || strcmp(quantity,'tcis') || strcmp(quantity,'vcis')
  
  % TODO: add oxygen
  
  if check_caa_sh_interval
    if ~exist('./.caa_sh_interval','file') && strcmp(quantity,'vcis')
      irf_log('proc','Inside magnetosphere. No vcis data fetched.')
      out_data = []; cd(old_pwd), return
    end
    if ~exist('./.caa_ms_interval','file') && strcmp(quantity,'edi')
      irf_log('proc','Outside magnetosphere. No edi data fetched.')
      out_data = []; cd(old_pwd), return
    end
  end
  
  ipref = '';
  r.qua = {quantity};
  switch(quantity)
    case 'b'
      r.ins = 'BPP';
      r.var = {'BPP'};
      r.qua = {'b'};
      
    case 'edi'
      r.ins = 'EDIR';
      r.var = {'EDI'};
      ipref = 'i'; % EDI pp is in the inertial frame
      
    case 'vcis'
      r.ins = 'CISR';
      r.qua = {'vcis_p', 'vcis_h'}; % CODIF and HIA
      r.var = {'VCp', 'VCh'};
      
    case 'ncis'
      r.ins = 'CISR';
      r.qua = {'ncis_p', 'ncis_h'}; % CODIF and HIA
      r.var = {'NCp', 'NCh'};
      
    case 'tcis'
      r.ins = 'CISR';
      r.qua = {'tcis_hpar','tcis_hper','tcis_ppar','tcis_pper'}; % CODIF and HIA
      r.var = {'TparCh', 'TperpCh','TparCp', 'TperpCp'};
      
    otherwise
      error('Check variable list')
  end
  
  save_file = ['./m' r.ins '.mat'];
  
  n_ok = 0;
  sax_loaded = 0;
  for i=1:length(r.qua)
    % first try ISDAT (fast) then files
    dat = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,r.qua{i});
    if isempty(dat) || ~any(any(~isnan(dat(:,2:end))))
      irf_log('dsrc',irf_ssub(...
        ['No data for ' ipref 'di' r.var{i} ', i' r.var{i} '?'],cl_id))
      continue
    end
    n_ok = n_ok + 1;
    c_eval([ipref r.var{i} '?=dat;'...
      'save_list=[save_list '' ' ipref r.var{i} '?''];'],cl_id);
    
    % Load SAX
    if ~sax_loaded
      sax_loaded = 1;
      [ok,sax] = c_load('SAX?',cl_id);
      if ~ok
        tempv = getData(cdb,start_time,dt,cl_id,'sax');
        if isempty(tempv)
          irf_log('dsrc',irf_ssub('Cannot load SAX?',cl_id))
          sax = [];
        else, sax = tempv{2};
        end
        clear tempv
      end
    end
    
    % Transform vector data to DSI
    if size(dat,2)>2
      if ~isempty(sax)
        c_eval([ipref 'di' r.var{i} '?=c_gse2dsi(dat,sax);'...
          'save_list=[save_list '' ' ipref 'di' r.var{i} '?''];'],cl_id);
      else
        irf_log('dsrc',...
          irf_ssub(['No data for ' ipref 'di' r.var{i} '?'],cl_id))
      end
    end
    clear dat
  end
  
  if ~n_ok, out_data = []; cd(old_pwd), return, end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % sax - spin axis orientation [GSE]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'sax')
  save_file='./mEPH.mat';
  
  % first try ISDAT (fast) then files
  lat = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slat');
  long = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slong');
  if ~isempty(lat), lat(isnan(lat(:,2)),:) = []; end
  if ~isempty(long), long(isnan(long(:,2)),:) = []; end
  if isempty(lat) || isempty(long)
    % Try ISDAT which does not handle data gaps
    [t,data] = caa_is_get(cdb.db, start_time, dt, ...
      cl_id, 'ephemeris', 'sax_lat');
    lat = [t data];
    [t,data] = caa_is_get(cdb.db, start_time, dt, ...
      cl_id, 'ephemeris', 'sax_long');
    long = [t data];
    
    if isempty(lat) || isempty(long)
      irf_log('dsrc',irf_ssub('No data for SAX?',cl_id))
      out_data = []; cd(old_pwd), return
    end
  end
  if (length(lat) ~= length(long)) || ~all(lat(:,1)==long(:,1))
    [ii1,ii2]=irf_find_comm_idx(lat,long);
    lat = lat(ii1,:); long = long(ii2,:);
  end
  if isempty(lat) || isempty(long)
    irf_log('dsrc',irf_ssub('No data for SAX?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  % XXX: TODO here we definitely need to check if the attitude is changing
  % with time instead of doing mean()
  lat = mean(lat(:,2)); long = mean(long(:,2));
  % Take first point only. This is OK according to AV
  [xspin,yspin,zspin] = sph2cart(long*pi/180,lat*pi/180,1);
  sax = [xspin yspin zspin]; %#ok<NASGU>
  
  eval(irf_ssub('SAX?=sax;save_list=[save_list '' SAX?''];',cl_id));
  clear sax lat long xspin yspin zspin
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % wbdwf - WBD waveforms.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'wbdwf')
  save_file = './mWBD.mat';
  try
    wf = c_wbd_read(start_time, dt, cl_id);  %#ok<NASGU>
    c_eval('wfWBD?=wf; save_list=[save_list '' wfWBD? ''];',cl_id);
  catch %#ok<CTCH>
    out_data = [];
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % whinat - WHISPER natural
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'whinat')
  save_file = './mWHI.mat';
  
  [t,data] = caa_is_get(cdb.db, start_time, dt, ...
    cl_id, 'whisper', 'natural');
  if isempty(data)
    irf_log('dsrc',irf_ssub('No data for WHINAT?',cl_id))
    out_data = []; cd(old_pwd), return
  end
  
  specrec.f = [1302.08	1464.84	1627.6	1790.36	1953.12	2115.89	2278.65	2441.41	2604.17	2766.93	2929.69	3092.45	3255.21	3417.97	3580.73	3743.49	3906.25	4069.01	4231.77	4394.53	4557.29	4720.05	4882.81	5045.57	5208.33	5371.09	5533.85	5696.61	5859.38	6022.14	6184.9	6347.66	6510.42	6673.18	6835.94	6998.7	7161.46	7324.22	7486.98	7649.74	7812.5	7975.26	8138.02	8300.78	8463.54	8626.3	8789.06	8951.82	9114.58	9277.34	9440.1	9602.86	9765.62	9928.39	10091.1	10253.9	10416.7	10579.4	10742.2	10904.9	11067.7	11230.5	11393.2	11556	11718.8	11881.5	12044.3	12207	12369.8	12532.6	12695.3	12858.1	13020.8	13183.6	13346.4	13509.1	13671.9	13834.6	13997.4	14160.2	14322.9	14485.7	14648.4	14811.2	14974	15136.7	15299.5	15462.2	15625	15787.8	15950.5	16113.3	16276	16438.8	16601.6	16764.3	16927.1	17089.8	17252.6	17415.4	17578.1	17740.9	17903.6	18066.4	18229.2	18391.9	18554.7	18717.4	18880.2	19043	19205.7	19368.5	19531.2	19694	19856.8	20019.5	20182.3	20345.1	20507.8	20670.6	20833.3	20996.1	21158.9	21321.6	21484.4	21647.1	21809.9	21972.7	22135.4	22298.2	22460.9	22623.7	22786.5	22949.2	23112	23274.7	23437.5	23600.3	23763	23925.8	24088.5	24251.3	24414.1	24576.8	24739.6	24902.3	25065.1	25227.9	25390.6	25553.4	25716.1	25878.9	26041.7	26204.4	26367.2	26529.9	26692.7	26855.5	27018.2	27181	27343.8	27506.5	27669.3	27832	27994.8	28157.6	28320.3	28483.1	28645.8	28808.6	28971.4	29134.1	29296.9	29459.6	29622.4	29785.2	29947.9	30110.7	30273.4	30436.2	30599	30761.7	30924.5	31087.2	31250	31412.8	31575.5	31738.3	31901	32063.8	32226.6	32389.3	32552.1	32714.8	32877.6	33040.4	33203.1	33365.9	33528.6	33691.4	33854.2	34016.9	34179.7	34342.4	34505.2	34668	34830.7	34993.5	35156.2	35319	35481.8	35644.5	35807.3	35970.1	36132.8	36295.6	36458.3	36621.1	36783.9	36946.6	37109.4	37272.1	37434.9	37597.7	37760.4	37923.2	38085.9	38248.7	38411.5	38574.2	38737	38899.7	39062.5	39225.3	39388	39550.8	39713.5	39876.3	40039.1	40201.8	40364.6	40527.3	40690.1	40852.9	41015.6	41178.4	41341.1	41503.9	41666.7	41829.4	41992.2	42154.9	42317.7	42480.5	42643.2	42806	42968.8	43131.5	43294.3	43457	43619.8	43782.6	43945.3	44108.1	44270.8	44433.6	44596.4	44759.1	44921.9	45084.6	45247.4	45410.2	45572.9	45735.7	45898.4	46061.2	46224	46386.7	46549.5	46712.2	46875	47037.8	47200.5	47363.3	47526	47688.8	47851.6	48014.3	48177.1	48339.8	48502.6	48665.4	48828.1	48990.9	49153.6	49316.4	49479.2	49641.9	49804.7	49967.4	50130.2	50293	50455.7	50618.5	50781.2	50944	51106.8	51269.5	51432.3	51595.1	51757.8	51920.6	52083.3	52246.1	52408.9	52571.6	52734.4	52897.1	53059.9	53222.7	53385.4	53548.2	53710.9	53873.7	54036.5	54199.2	54362	54524.7	54687.5	54850.3	55013	55175.8	55338.5	55501.3	55664.1	55826.8	55989.6	56152.3	56315.1	56477.9	56640.6	56803.4	56966.1	57128.9	57291.7	57454.4	57617.2	57779.9	57942.7	58105.5	58268.2	58431	58593.8	58756.5	58919.3	59082	59244.8	59407.6	59570.3	59733.1	59895.8	60058.6	60221.4	60384.1	60546.9	60709.6	60872.4	61035.2	61197.9	61360.7	61523.4	61686.2	61849	62011.7	62174.5	62337.2	62500	62662.8	62825.5	62988.3	63151	63313.8	63476.6	63639.3	63802.1	63964.8	64127.6	64290.4	64453.1	64615.9	64778.6	64941.4	65104.2	65266.9	65429.7	65592.4	65755.2	65918	66080.7	66243.5	66406.2	66569	66731.8	66894.5	67057.3	67220.1	67382.8	67545.6	67708.3	67871.1	68033.9	68196.6	68359.4	68522.1	68684.9	68847.7	69010.4	69173.2	69335.9	69498.7	69661.5	69824.2	69987	70149.7	70312.5	70475.3	70638	70800.8	70963.5	71126.3	71289.1	71451.8	71614.6	71777.3	71940.1	72102.9	72265.6	72428.4	72591.1	72753.9	72916.7	73079.4	73242.2	73404.9	73567.7	73730.5	73893.2	74056	74218.8	74381.5	74544.3	74707	74869.8	75032.6	75195.3	75358.1	75520.8	75683.6	75846.4	76009.1	76171.9	76334.6	76497.4	76660.2	76822.9	76985.7	77148.4	77311.2	77474	77636.7	77799.5	77962.2	78125	78287.8	78450.5	78613.3	78776	78938.8	79101.6	79264.3	79427.1	79589.8	79752.6	79915.4	80078.1	80240.9	80403.6	80566.4	80729.2	80891.9	81054.7	81217.4	81380.2	81543];
  specrec.p = {data'};
  specrec.t = t; %#ok<STRNU>
  
  c_eval('WHINAT?=specrec;save_list=[save_list '' WHINAT? ''];',cl_id);
  clear t data specrec
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else, error('caa:noSuchQuantity','Quantity ''%s'' unknown',quantity)
end %main QUANTITY

% saving
% If flag_save is set, save variables to specified file
if flag_save==1 && ~isempty(save_list) && ~isempty(save_file)
  irf_log('save',[save_list ' -> ' save_file])
  if exist(save_file,'file')
    eval(['save -append ' save_file ' ' save_list]);
  else
    eval(['save ' save_file ' ' save_list]);
  end
end

% prepare the output
if nargout > 0
  if ~isempty(save_list)
    sl = tokenize(save_list);
    out_data = {sl};
    for i=1:length(sl)
      eval(['out_data{i+1}=' sl{i} ';'])
    end
  end
else
  clear out_data
end

cd(old_pwd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=check_timeline(data)
% This function is aimed at finding negative time jumps
% 1) if it is a single misplaced pachet - correct it
% 2) otherwise blank +/- 10 sec of data around the jump
% ISDAT BUG 11 http://squid.irfu.se/bugzilla/show_bug.cgi?id=11

out = data;

% Find jumps back larger than BAD_DT sec.
BAD_DT = 0;
% Remove BAD_MARGIN from each side
BAD_MARGIN = 4.5;
MAX_JITTER = 0.0005;

while 1
  iNegative = find(diff(out(:,1))<=-BAD_DT);
  if isempty(iNegative), return, end
  iNegative = iNegative(1);
  
  tbj = out(iNegative,1) + BAD_MARGIN;
  taj = out(iNegative+1,1) - BAD_MARGIN;
  idx = find(out(:,1) > taj & out(:,1) < tbj);
  ttmp = out(idx,1); ddt = diff(ttmp); mddt = median(ddt);
  iiPositive = find( ddt - mddt > mddt*MAX_JITTER );
  
  i1secJump = [];
  if ~isempty(iiPositive)
    dtJump = abs(out(iNegative+1,1)-ttmp(iiPositive));
    nSec = round(dtJump);
    idx1secOK = abs(dtJump-nSec)<0.05;
    i1secJump = iiPositive(idx1secOK) + idx(1) -1; % absolute index
    if any(idx1secOK)
      if length(i1secJump) > 1
        nSec = nSec(idx1secOK);
        i1secJump = i1secJump(nSec == min(nSec));
        i1secJump = i1secJump(1); nSec = min(nSec);
      end
      
      dtNeg = out(iNegative+1,1)-out(iNegative,1);
      dtPos = out(i1secJump+1,1) - out(i1secJump,1);
      if abs(dtNeg)>1e-6 && abs(dtNeg+dtPos)>3*mddt % expected to diff by 2 nominal time steps
        i1secJump = [];
      else
        if abs(dtNeg)<1e-6, dt = -mddt; % shift by one data point
        else, dt = (dtNeg-dtPos)/2;
        end
        if iNegative>i1secJump, iCorr = (i1secJump+1):iNegative; % shift backward
        else, dt = -dt; iCorr = (iNegative+1):i1secJump; % shift forward
        end
        out(iCorr,1) = out(iCorr,1) + dt;
        irf_log('proc',...
          sprintf('Correcting %d misplaced packet(s) at %s by %.4f sec',...
          nSec,epoch2iso(out(iCorr(1),1)),dt));
      end
    end
  end
  
  if isempty(i1secJump)
    out(out(:,1) > taj & out(:,1) < tbj,:) = [];
    irf_log('proc',...
      sprintf('Bad time %s - %s',epoch2iso(taj,1),epoch2iso(tbj,1)));
  end
  
  if isempty(out), return, end % extra precaution to avoid a dead loop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = prob_s(ns_ops_rec,warn)
if nargin<2, warn=0; end
if warn, s = 'WARNING: ';
else, s = 'PROBLEM: ';
end
ss = [s caa_errid2str(ns_ops_rec(4)) ' ' epoch2iso(ns_ops_rec(1),1)...
  ' -- ' epoch2iso(ns_ops_rec(1)+ns_ops_rec(2),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataout = rm_ib_spike(data)
% Remove spikes in staff IB data

DT = 0.455; % sec, time window
THR = 5; % data above TRH StDev discarded

dataout = data;

ndata = length(data(:,1));
dt2 = ceil(DT*c_efw_fsample(data,'ib')/2);

i = 1;
while i < ndata
  i2 = i + dt2*2;
  if i2 > ndata % last window
    i = ndata - dt2*2;
    i2 = ndata;
  end
  x = data(i:i2,2:4);
  y = x;
  
  x = detrend(x);
  s = std(x);
  for comp=3:-1:1
    ii = find(abs(x(:,comp))>THR*s(comp));
    if ~isempty(ii), y(ii,:) = y(ii+1,:); end
  end
  
  dataout(i:i2,2:4) = y;
  
  i = i2 + 1;
end

