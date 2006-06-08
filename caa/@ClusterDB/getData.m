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
%	pburst: P{4kHz,32kHz}{cl_id}p{1..4}, wbE{cl_id}p12,34 -> mEFWburst	
%			// probe potentials (4kHz,32kHz), and electric fields
%
%	//// Ephemeris ////
%	sax : SAX{cl_id} ->mEPH
%			// spin axis vector [GSE] 
%	a   : A{cl_id} -> mA	// SC phase
%	r   : R{cl_id} -> mR	// SC position
%	v   : V{cl_id}, diV{cl_id} -> mR	// SC velocity
%
%	//// Other instruments ////
%	b   : BPP{cl_id},diBPP{cl_id}	->mBPP	// B FGM PP [GSE+DSI] 
%	bfgm: B{cl_id},diB{cl_id}	->mB	// B FGM** [GSE+DSI]
%		** contact Stephan Buchert
%	edi : iEDI{cl_id},idiEDI{cl_id}	->mEDI	// E EDI PP (inert frame) [GSE+DSI] 
%	ncis: NC(p,h){cl_id}			->mCIS	// N CIS PP
%	tcis: T(par,perp)C(p,h){cl_id}	->mCIS	// T CIS PP 
%	vcis: VC(p,h){cl_id},diVC(p,h){cl_id}  ->mCIS	// V CIS PP [GSE+DSI] 
%	wbdwf: wfWBD{cl_id} -> mWBD	// WBD waveforms E/B 
%
%	options - one of the following:
%	nosave : do no save on disk
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
%
% $Id$

% Copyright 2004,2005 Yuri Khotyaintsev
% Parts of the code are (c) Andris Vaivads

error(nargchk(5,15,nargin))
if nargin > 5, property_argin = varargin; end

out_data = '';

% default options
flag_save = 1;

for i=1:length(varargin)
	switch(varargin{i})
	case 'nosave'
		flag_save = 0;
	otherwise
		irf_log('fcal',['Option ''' varargin{i} '''not recognized'])
	end
end

save_file = '';
save_list = '';

old_pwd = pwd;

%Create the storage directory if it does not exist
if ~exist(cdb.sp, 'dir')
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir(cdb.sp);
	if SUCCESS, irf_log('save',['Created storage directory ' cdb.sp])
    else error(MESSAGE)
	end
end

cd(cdb.sp) %enter the storage directory
irf_log('save',['Storage directory is ' cdb.sp])

% Create .interval
if ~exist('./.interval','file')
	fid = fopen('.interval','w');
	if fid<0, irf_log('save','problem creating .interval'),cd(old_pwd),return, end
	count = fprintf(fid,'%s %s',epoch2iso(start_time),num2str(dt));	fclose(fid);
	if count<=0, irf_log('save','problem writing to .interval'),cd(old_pwd),return, end
end

% Read list of nonstandard operations and see if we have one of those 
% during the requested period. Permanent problems (as loss of 
% probes, filters, etc.) must be programmed separately
if strcmp(quantity,'e') || strcmp(quantity,'eburst') ||...
	strcmp(quantity,'p') || strcmp(quantity,'pburst')
	
	ns_ops = c_ctl('get',cl_id,'ns_ops');
	if isempty(ns_ops)
		c_ctl('load_ns_ops',[cdb.dp '/caa-control'])
		ns_ops = c_ctl('get',cl_id,'ns_ops');
	end
	if ~isempty(ns_ops)
		if strcmp(quantity,'p') || strcmp(quantity,'pburst')
			errlist = caa_str2errid('hxonly');
		else errlist = [];
		end
		[start_time_nsops, dt_nsops] = caa_ns_ops_int(start_time,dt,ns_ops,errlist);
		if isempty(start_time_nsops)
			irf_log('dsrc',sprintf('bad NS_OPS interval for C%d - %s',cl_id,quantity))
			data = []; 
			cd(old_pwd), return
		end
    else start_time_nsops = start_time; dt_nsops = dt;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dsc - EFW DSC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'dsc')
	save_file = './mEFWR.mat';
	
	[t,dsc] = caa_is_get(cdb.db, start_time, dt, cl_id, 'efw', 'DSC');
	if isempty(dsc)
		irf_log('dsrc',irf_ssub('No data for DSC?',cl_id))
		data = []; cd(old_pwd), return
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
				t_start_save(l+1) = t(i);
				dsc_save(:,l+1) = dsc(dsc_is,i);
				jump_flag(l+1) = 1;
				n_jumpy = n_jumpy + 1;
			else
				% Save a good interval
				t_start_save(l+1) = t_st;
				dsc_save(:,l+1) = dsc_good(dsc_is);
				jump_flag(l+1) = 0;
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
		t_start_save(end+1) = t_st;
		dsc_save(:,end+1) = dsc_good(dsc_is);
		n_good = n_good + 1;					
		irf_log('dsrc',['Saving good from ' epoch2iso(t_st,1) '-' epoch2iso(t_end,1)])
	end
	
	irf_log('dsrc',sprintf('\nFound total %d good and %d jumpy intervals',...
		n_good, n_jumpy))
	c_eval('DSC?=[t_start_save dsc_save''];save_list=[save_list '' DSC? ''];',cl_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fdm - EFW FDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'fdm')
	save_file = './mEFWR.mat';
	
	[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, 'efw', 'FDM');
	if isempty(data)
		irf_log('dsrc',irf_ssub('No data for FDM?',cl_id))
		data = []; cd(old_pwd), return
    else c_eval('FDM?=[t data''];',cl_id);
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
	
	% Check for p1 problems on SC1 and SC3
	if (start_time>toepoch([2001 12 28 03 00 00]) && cl_id==1) || ...
		(start_time>toepoch([2002 07 29 09 06 59 ]) && cl_id==3)
		probe_list = 2:4;
		irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
	end
	
	for probe=probe_list;
		[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, ...
			'efw', 'E', ['p' num2str(probe)],'bias');
		if isempty(data)
			irf_log('dsrc',irf_ssub('No data for IBIAS?p!',cl_id,probe))
		else
			eval(irf_ssub('IBIAS?p!=[t data];',cl_id,probe))
			p_ok = [p_ok probe];
		end
		clear t data
	end
	
	if isempty(p_ok), data = []; cd(old_pwd), return, end
	for probe=p_ok;
		eval(irf_ssub('save_list=[save_list ''IBIAS?p! ''];',cl_id,probe)) 
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% efwt - EFW clock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'efwt')
	save_file = './mEFWR.mat';
	
	% Read EFW clock to check for time since last reset
	t = []; data = []; efwtime = [];
	for st_tmp = start_time-16:32:start_time+dt+16
		[t_tmp,data] = caa_is_get(cdb.db, st_tmp, 32, cl_id, 'efw', 'DSC');
		if ~isempty(data)
			efwtime_tmp = (data(81,:) +data(82,:)*256 +data(83,:)*65536 + ...
				data(84,:)*16777216 +data(85,:)*4294967296)/1000;
			efwtime = [efwtime efwtime_tmp];
			t = [t t_tmp'];
		end
	end
	
	if isempty(efwtime)
		irf_log('dsrc',irf_ssub('No data for EFWT?',cl_id))
		data = []; cd(old_pwd), return
	end
	c_eval(['EFWT?=[t; efwtime]'';'...
		'save_list=[save_list '' EFWT? ''];'],cl_id);
	clear t data efwtime
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmode - EFW tape mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'tmode')
	save_file = './mEFWR.mat';
	
	%% Find TapeMode
	% We read FDM from isdat and 5-th column contains the HX mode
	% (undocumented feature)
	% 0 - normal mode  (V12L,V34L)
	% 1 - tape mode 1  (V12M,V34M)
	% 2 - tape mode 2  (V12M,V34M)
	% 3 - tape mode 3  (V1M,V2M,V3M,V4M)
	[t,data] = caa_is_get(cdb.db,start_time,dt,cl_id,'efw','FDM');
	if isempty(data)
		irf_log('dsrc',irf_ssub('No data for mTMode?',cl_id))
		data = []; cd(old_pwd), return
	end
	
	c_eval(['mTMode?=data(5,:);'...
		'save_list=[save_list '' mTMode? ''];'],cl_id)
	clear t data
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e - Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'e') || strcmp(quantity,'eburst')
	
	if strcmp(quantity,'eburst'), do_burst = 1; else do_burst = 0; end
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
				data = []; cd(old_pwd), return
			end
			tm = tmode{2};
			clear tmode
		end
		if isempty(tm)
			irf_log('dsrc',irf_ssub('Cannot load mTMode?',cl_id))
			data = []; cd(old_pwd), return
		end
			
		if any(tm~=tm(1))
			irf_log('dsrc','tape mode changes during the selected time inteval')
			irf_log('dsrc','data interval will be truncated')
		end
		tm = tm(1);
		if tm<1e-30, param='10Hz'; else param='180Hz'; end
		clear tm
		
		%%%%%%%%%%%%%%%%%%%%%%%%% FILTER MAGIC %%%%%%%%%%%%%%%%%%%%%%
		if (((cl_id==1 && start_time>toepoch([2001 07 30 17 05 54.9])) || ...
			(cl_id==3 && start_time>toepoch([2001 07 31 00 12 29.5]))) && ...
			start_time<toepoch([2001 09 02 23 15 00])) || ...
			(cl_id==4 && ((start_time>toepoch([2001 07 31 04 55 33.15]) && ...
			start_time<toepoch([2001 08 02 11 25 40])) || ...
			(start_time>toepoch([2001 08 06 23 58 50.7]) && ...
			start_time<toepoch([2001 09 02 23 15 00]))))
			% all sc run on 180Hz filter in august 2001 most of the time
			param='180Hz';
		elseif start_time>toepoch([2001 09 10 04 21 57.6]) && ...
			start_time<toepoch([2001 09 15 06 30 00])
			% this needs to be investigated.... 
			param='180Hz';
		elseif cl_id==2 && start_time>toepoch([2001 07 23 00 00 00])
			% 10Hz filter problem on SC2
			param='180Hz';
		end
		%%%%%%%%%%%%%%%%%%%%%%% END FILTER MAGIC %%%%%%%%%%%%%%%%%%%%
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%% PROBE MAGIC %%%%%%%%%%%%%%%%%%%%%%
	pl = [12,34];
	if (cl_id==1 || cl_id==3) && (start_time>toepoch([2003 9 29 00 27 0]) || ...
		(start_time>toepoch([2003 3 27 03 50 0]) && start_time<toepoch([2003 3 28 04 55 0])) ||...
		(start_time>toepoch([2003 4 08 01 25 0]) && start_time<toepoch([2003 4 09 02 25 0])) ||...
		(start_time>toepoch([2003 5 25 15 25 0]) && start_time<toepoch([2003 6 08 22 10 0])) )
		
		% FSW 2.4. Use P32 on SC1 and SC3
		pl = [32, 34];
		irf_log('dsrc',sprintf('            !Using p32 for sc%d',cl_id));
		
	elseif (cl_id==1 && start_time>toepoch([2001 12 28 03 00 00])) || ...
		(cl_id==3 && start_time>toepoch([2002 07 29 09 06 59 ]))
		
		% p1 problems on SC1 and SC3
		pl = 34;
		irf_log('dsrc',sprintf('            !Only p34 exists for sc%d',cl_id));
	end
	%%%%%%%%%%%%%%%%%%%%%%% END PROBE MAGIC %%%%%%%%%%%%%%%%%%%%
	
	p_ok = [];
	for probe=pl
		irf_log('dsrc',['EFW...sc' num2str(cl_id)...
			'...Ep' num2str(probe) ' ' param ' filter']);
		t = [];	data = [];
		for in=1:length(start_time_nsops)
			if length(start_time_nsops)>1
				irf_log('dsrc',...
					sprintf('chunk #%d : %s %d sec',in,...
						epoch2iso(start_time_nsops(in),1),dt_nsops(in)))
			end
			[t_tmp,data_tmp] = caa_is_get(cdb.db, start_time_nsops(in), dt_nsops(in), cl_id, ...
				'efw', 'E', ['p' num2str(probe)], param, tmmode);
			t = [t; t_tmp]; data = [data; data_tmp]; clear t_tmp data_tmp
		end
		if ~isempty(data)
			% Correct start time of the burst
			if do_burst
				burst_f_name = irf_ssub([irf_fname(t(1),1) 'we.0?'],cl_id);
				burst_f_name = [cdb.dp '/burst/' burst_f_name];
				if exist(burst_f_name,'file')
					db = Mat_DbOpen(cdb.db);
					err_t = t(1) - c_efw_burst_chkt(db,burst_f_name);
					Mat_DbClose(db);
					
					t = t - err_t;
					irf_log('dsrc',...
						['burst start time was corrected by ' num2str(err_t) ' sec'])
                else irf_log('dsrc','burst start time was not corrected')
				end
			end
			
			data = check_timeline([t data]);
			c_eval([var_name  num2str(probe) '=data;'],cl_id);
			p_ok = [p_ok probe];
			
		else
			irf_log('dsrc',...
				irf_ssub(['No data for ' var_name num2str(probe)],cl_id))
		end
	end
	
	if isempty(p_ok), data = []; cd(old_pwd), return, end
	for probe=p_ok
		c_eval(['save_list=[save_list '' ' var_name num2str(probe) '''];'],cl_id);
	end
	clear t data tm pl param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p - SC potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'p') || strcmp(quantity,'pburst')
	
	if strcmp(quantity,'pburst'), do_burst = 1; else do_burst = 0; end
	if do_burst 
		save_file = './mEFWburstR.mat';
		tmmode='burst';
		param={'4kHz','32kHz'};
		var_name = 'wbE?p';
	else
		save_file = './mPR.mat';
		param={'10Hz'}; tmmode='lx';
	end
	
	probe_list = 1:4;
	
	%%%%%%%%%%%%%%%%%%%%%%%%% PROBE MAGIC %%%%%%%%%%%%%%%%%%%%%%
	% Check for p1 problems on SC1 and SC3
	if (cl_id==1 && start_time>toepoch([2001 12 28 03 00 00])) || ...
		(cl_id==3 && start_time>toepoch([2002 07 29 09 06 59 ]))
		probe_list = 2:4;
		p1 = [];
		irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
	end
	% 10Hz filter problem on C2 p3
	% Any changes should also go to ClusterProc/getData/probesa
	if cl_id==2 && start_time>toepoch([2001 07 23 00 00 00]) && ~do_burst
		probe_list = [1 2 4];
		p3 = [];
		irf_log('dsrc',sprintf('10Hz filter problem on sc%d',cl_id));
	end
	%%%%%%%%%%%%%%%%%%%%%%% END PROBE MAGIC %%%%%%%%%%%%%%%%%%%%
	
	n_ok = 0;
	for j=1:length(param), for probe=probe_list;
    	irf_log('dsrc',['EFW...sc' num2str(cl_id) '...probe' num2str(probe)...
			'->P' param{j} num2str(cl_id) 'p' num2str(probe)]);
		t = [];	data = [];
		for in=1:length(start_time_nsops)
			if length(start_time_nsops)>1
				irf_log('dsrc',...
					sprintf('chunk #%d : %s %d sec',in,...
						epoch2iso(start_time_nsops(in),1),dt_nsops(in)))
			end
			[t_tmp,data_tmp] = caa_is_get(cdb.db, start_time_nsops(in), dt_nsops(in), cl_id, ...
				'efw', 'E', ['p' num2str(probe)],param{j}, tmmode);
			t = [t; t_tmp]; data = [data; data_tmp]; clear t_tmp data_tmp
		end
		
		if ~isempty(data)
			% Correct start time of the burst
			if do_burst
				burst_f_name = irf_ssub([irf_fname(t(1),1) 'we.0?'],cl_id);
				burst_f_name = [cdb.dp '/burst/' burst_f_name];
				if exist(burst_f_name,'file')
					db = Mat_DbOpen(cdb.db);
					err_t = t(1) - c_efw_burst_chkt(db,burst_f_name);
					irf_log('dsrc',['burst start time was corrected by ' ...
						num2str(err_t) ' sec'])
					Mat_DbClose(db);
					t = t - err_t;
				else
					irf_log('dsrc','burst start time was not corrected')
				end
			end
			
			data = check_timeline([t data]);
			eval(irf_ssub(['P' param{j} '?p!=data;'...
				'save_list=[save_list ''P' param{j} '?p! ''];'],cl_id,probe));
			n_ok = n_ok + 1;
			
        else irf_log('dsrc', irf_ssub(['No data for P' param{j} '?p!'],cl_id,probe));
		end 
		clear t data
    end, end
	
    if ~n_ok, data = []; cd(old_pwd), return, end
	
	% Make electric field for the burst
	% TODO: move to ClusterProc
	if do_burst
		for j=1:length(param)
			for probe=[1 3]
				if exist(irf_ssub(['P' param{j} '?p!'],cl_id,probe),'var') && ...
				exist(irf_ssub(['P' param{j} '?p!'],cl_id,probe+1),'var')
					eval(irf_ssub(['E(:,1)=P' param{j} '?p$(:,1);E(:,2)=(P' ...
						param{j} '?p$(:,2)-P' param{j} '?p!(:,2))/.088;'],...
						cl_id,probe,probe+1));
					vn = [var_name num2str(probe) num2str(probe+1)];
					if exist(irf_ssub(vn,cl_id),'var')
						c_eval(['tmpE=' vn ';'],cl_id)
						if tmpE(1,1) > E(1,1)
							tmpE(:,end+1:end+size(E,1)) = E;
							E = tmpE;
						else
							E(:,end+1:end+size(tmpE,1)) = tmpE;
						end
						clear tmpE
					end
					c_eval([vn '=E;save_list=[save_list ''' vn  ' ''];'],cl_id);
					clear E
				end
			end
		end
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Phase, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'a')
	save_file = './mA.mat';
	
	n_ok = 0;
	
	% We ask for 2 sec more from each side 
	% to avoid problems with interpolation.
	[t,data] = caa_is_get(cdb.db, start_time-2, dt+4, ...
		cl_id, 'ephemeris', 'phase');
	if ~isempty(data)
		c_eval('A?=[t data];save_list=[save_list '' A? ''];',cl_id);
		n_ok = n_ok + 1;
    else irf_log('dsrc',irf_ssub('No data for A?',cl_id))
	end
	clear t data
	
	[t,data] = caa_is_get(cdb.db, start_time-2, dt+4, ...
		cl_id, 'ephemeris', 'phase_2');
	if ~isempty(data)
		c_eval('Atwo?=[t data];save_list=[save_list '' Atwo? ''];',cl_id);
		n_ok = n_ok + 1;
	else
		c_eval('Atwo?=[];',cl_id);
		irf_log('dsrc',irf_ssub('No data for Atwo?',cl_id))
	end
	clear t data
	
	if ~n_ok, data = []; cd(old_pwd), return, end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'r')
	save_file = './mR.mat';
	
	[t,data] = caa_is_get(cdb.db, start_time, dt, ...
		cl_id, 'ephemeris', 'position');
	if ~isempty(data), c_eval('R?=[t data''];save_list=[save_list '' R? ''];',cl_id);
	else
		irf_log('dsrc',irf_ssub('No data for R?',cl_id))
		data = []; cd(old_pwd), return
	end
	clear t data
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'v')
	save_file = './mR.mat';
	
	[t,data] = caa_is_get(cdb.db, start_time, dt, ...
		cl_id, 'ephemeris', 'velocity');
	if isempty(data)
		irf_log('dsrc',irf_ssub('No data for V?, diV?',cl_id))
		data = []; cd(old_pwd), return
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
        else sax = tempv{2};
		end
		clear tempv
	end
	if ~isempty(sax)
		c_eval('diV?=c_gse2dsi(V?,sax);save_list=[save_list '' diV? ''];',cl_id);
    else irf_log('dsrc',irf_ssub('No data for diV?',cl_id))
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
% B FGM - full res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'bfgm')
	save_file = './mB.mat';
	
	dat = c_get_bfgm(start_time + [0 dt],cl_id);
	
	if isempty(dat)
		irf_log('dsrc',irf_ssub('No data for B, diB?',cl_id))
		data = []; cd(old_pwd), return
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
        else sax = tempv{2};
		end
		clear tempv
	end
	if ~isempty(sax)
		c_eval('diB?=c_gse2dsi(B?,sax);save_list=[save_list '' diB? ''];',cl_id);
    else irf_log('dsrc',irf_ssub('No data for diB?',cl_id))
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CSDS PP [GSE+DSI] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'b') || strcmp(quantity,'edi') || ...
	strcmp(quantity,'ncis') || strcmp(quantity,'tcis') || strcmp(quantity,'vcis')
	
	% TODO: add oxygen
	
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
		if isempty(dat)
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
                else sax = tempv{2};
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
	
	if ~n_ok, data = []; cd(old_pwd), return, end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sax - spin axis orientation [GSE] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'sax')
	save_file='./mEPH.mat';

	% first try ISDAT (fast) then files
	lat = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slat');
	long = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slong');
	
	if isempty(lat) || isempty(long)
		irf_log('dsrc',irf_ssub('No data for SAX?',cl_id))
		data = []; cd(old_pwd), return
	end
	
	% Take first point only. This is OK according to AV
	[xspin,yspin,zspin] = sph2cart(long(1,2)*pi/180,lat(1,2)*pi/180,1);
	sax = [xspin yspin zspin];

	eval(irf_ssub('SAX?=sax;save_list=[save_list '' SAX?''];',cl_id));
	clear sax lat long xspin yspin zspin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wbdwf - WBD waveforms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'wbdwf')
	save_file = './mWBD.mat';
	try
		wf = c_wbd_read(start_time, dt, cl_id);
		c_eval('wfWBD?=wf;',cl_id); 
		c_eval('save_list=[save_list '' wfWBD? ''];',cl_id);
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else, error('caa:noSuchQuantity','Quantity ''%s'' unknown',quantity)
end %main QUANTITY

% saving
% If flag_save is set, save variables to specified file
if flag_save==1 && ~isempty(save_list) && length(save_file)>0
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
% This function is aimed at removing pieces of data
% with time jumping back
% ISDAT BUG 11 http://squid.irfu.se/bugzilla/show_bug.cgi?id=11

out = data;

% Find jumps back larger than BAD_DT sec.
BAD_DT = 0.5;
% Remove BAD_MARGIN from each side
BAD_MARGIN = 10;

while 1
	ii = find(diff(out(:,1))<-BAD_DT);
	if isempty(ii), return, end
	
	tbj = out(ii(1),1) + BAD_MARGIN;
	taj = out(ii(1)+1,1) - BAD_MARGIN;	
	out(out(:,1) > taj & out(:,1) < tbj,:) = [];
	irf_log('proc',...
		sprintf('Bad time %s - %s',epoch2iso(taj,1),epoch2iso(tbj,1)));
		
	if isempty(out), return, end % extra precaution to avoid a dead loop
end

