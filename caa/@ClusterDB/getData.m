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
%	e   : wE{cl_id}p12,34 -> mER
%			// electric fields (HX)
%	p   : P{cl_id} -> mPR	
%			// probe potential (LX)
%	dsc	: DSC{cl_id} -> mFDM
%			// EFW DSC
%	fdm	: FDM{cl_id} -> mFDM
%			// EFW FDM
%	ibias: IBIAS{cl_id}p{1..4} -> mFDM
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

start_date_str = strrep(datestr(fromepoch(start_time),29),'-','');

save_file = '';
save_list = '';

old_pwd = pwd;

%Create the storage directory if it does not exist
if ~exist(cdb.sp, 'dir')
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir(cdb.sp);
	if SUCCESS, irf_log('save',['Created storage directory ' cdb.sp])
	else, error(MESSAGE)
	end
end

cd(cdb.sp) %enter the storage directory
irf_log('save',['Storage directory is ' cdb.sp])

% Create .interval
if ~exist('./.interval','file')
	[s,w] = unix(...
		['echo "' epoch2iso(start_time) ' ' num2str(dt) '">.interval'],'-echo');
	if s~=0, irf_log('save','problem creating .interval'),cd(old_pwd),return, end
	irf_log('save','creating .interval')
end

% Read list of nonstandard operations and see if we have one of those 
% during the requested period
if strcmp(quantity,'e')|strcmp(quantity,'eburst')|strcmp(quantity,'p')|strcmp(quantity,'pburst')
	ns_ops = c_ctl('get',cl_id,'ns_ops');
	if isempty(ns_ops)
		c_ctl('load_ns_ops',[cdb.dp '/caa'])
		ns_ops = c_ctl('get',cl_id,'ns_ops');
	end
	if ~isempty(ns_ops)
		% remove records which cover permanent problems (as loss of 
		% probes, filters, etc.) as these must be programmed separately
		ii = find(ns_ops(:,2)==-1);
		ns_ops(ii,:) = [];
		
		ii = find(ns_ops(:,1)<start_time+dt & ns_ops(:,1)>start_time);
		if ~isempty(ii)
			for j=1:length(ii)
				if ns_ops(ii(j),4)<10
					% no/bad data - truncate the period
					irf_log('proc',...
					['PROBLEM: ' caa_errid2str(ns_ops(ii(j),4)) ' from ' ...
					epoch2iso(ns_ops(ii(j),1)) ' to ' ...
					epoch2iso(ns_ops(ii(j),1)+ns_ops(ii(j),2))])
					irf_log('proc',...
					['setting dt to ' num2str(ns_ops(ii(j),1)-start_time)])
					dt = ns_ops(ii(j),1) - start_time;
				else
					irf_log('proc',...
					['WARNING: ' caa_errid2str(ns_ops(ii(j),4)) ' from ' ...
					epoch2iso(ns_ops(ii(j),1)) ' to ' ...
					epoch2iso(ns_ops(ii(j),1)+ns_ops(ii(j),2))])
				end
			end
		end
		% clear already processed records
		ns_ops(ii,:) = [];
		ii = find(ns_ops(:,1)+ns_ops(:,2)<start_time+dt & ns_ops(:,1)+ns_ops(:,2)>start_time);
		if ~isempty(ii)
			for j=1:length(ii)
				if ns_ops(ii(j),4)<10
					% no/bad data - truncate the period
					irf_log('proc',...
					['PROBLEM: ' caa_errid2str(ns_ops(ii(j),4)) ' from ' ...
					epoch2iso(ns_ops(ii(j),1)) ' to ' ...
					epoch2iso(ns_ops(ii(j),1)+ns_ops(ii(j),2))])
					irf_log('proc',...
					['setting start_time to ' num2str(ns_ops(ii(j),1))])
					start_time = ns_ops(ii(j),1);
				else
					irf_log('proc',...
					['WARNING: ' caa_errid2str(ns_ops(ii(j),4)) ' from ' ...
					epoch2iso(ns_ops(ii(j),1)) ' to ' ...
					epoch2iso(ns_ops(ii(j),1)+ns_ops(ii(j),2))])
				end
			end
		end
	end
end

% We need to chech whether FDM was fetched, and wherether SWEEP and BDUMP 
% are computed
if strcmp(quantity,'e')|strcmp(quantity,'p')
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dsc - EFW DSC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'dsc')
	save_file = './mFDM.mat';
	
	[t,dsc] = caa_is_get(cdb.db, start_time, dt, cl_id, 'efw', 'DSC');
	if isempty(dsc)
		irf_log('dsrc',irf_ssub('No data for DSC?',cl_id))
		return
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
	t_stop_save = [];
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
			% the same as previous
			count_good = count_good + 1;
			if count_good == 1
				t_st = t_dsc_last;
				dsc_good = dsc_last;
				%irf_log('dsrc',['Found good at ' epoch2iso(t_dsc_last,1)])
			end
			t_end = t(i);
		else
			% differ from previous
			dsc_last = dsc(:,i);
			l = length(jump_flag);
			if count_good == 0
				%the previous point was also different
				t_start_save(l+1) = t(i);
				%t_stop_save(l+1) = t(i);
				dsc_save(:,l+1) = dsc(dsc_is,i);
				jump_flag(l+1) = 1;
				n_jumpy = n_jumpy + 1;
			else
				% save the good interval
				t_start_save(l+1) = t_st;
				%t_stop_save(l+1) = t_end;
				dsc_save(:,l+1) = dsc_good(dsc_is);
				jump_flag(l+1) = 0;
				n_good = n_good + 1;
				irf_log('dsrc',['Saving good from ' epoch2iso(t_st,1) '-' epoch2iso(t_end,1)])
			end
			count_good = 0;
		end
		t_dsc_last = t(i);
	end

	if count_good > 0 %last interval was also good
		% save the good interval
		t_start_save(end+1) = t_st;
		%t_stop_save(end+1) = t_end;
		dsc_save(:,end+1) = dsc_good(dsc_is);
		n_good = n_good + 1;					
		irf_log('dsrc',['Saving good from ' epoch2iso(t_st,1) '-' epoch2iso(t_end,1)])
	end
	
	irf_log('dsrc',sprintf('\nFound total %d good and %d jumpy intervals',...
		n_good, n_jumpy))
	c_eval(['DSC?=[t_start_save dsc_save''];save_list=[save_list '' DSC? ''];'],cl_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fdm - EFW FDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'fdm')
	save_file = './mFDM.mat';
	
	[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, 'efw', 'FDM');
	if isempty(data)
		irf_log('dsrc',irf_ssub('No data for FDM?',cl_id))
		return
	end
	c_eval(['FDM?=[t data''];save_list=[save_list '' FDM? ''];'],cl_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ibias - EFW probe bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'ibias')
	save_file = './mFDM.mat';
	
	probe_list = 1:4;
	
	% Check for p1 problems on SC1 and SC3
	if (start_time>toepoch([2001 12 28 03 00 00])&cl_id==1) | (start_time>toepoch([2002 07 29 09 06 59 ])&cl_id==3)
		probe_list = 2:4;
		p1 = [];
		irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
	end
	
	for probe=probe_list;
		[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, ...
		'efw', 'E', ['p' num2str(probe)],'bias');
		if ~isempty(data)
			data = double(real(data));
			t = double(t);
			eval(irf_ssub(...
				['p!=[t data];save_list=[save_list ''IBIAS?p! ''];IBIAS?p!=p!;'],...
				cl_id,probe)) 
			clear t data
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e - Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'e') | strcmp(quantity,'eburst')
	
	if strcmp(quantity,'eburst'), do_burst = 1; else do_burst = 0; end
	if do_burst 
		save_file = './mEFWburst.mat';
		tmmode='burst';
		param='8kHz';
		var_name = 'wbE?p';
	else 
		save_file = './mER.mat';
		tmmode='hx';
		var_name = 'wE?p';
	
		%% Find TapeMode
		% We read FDM from isdat and 5-th column contains the HX mode
		% (undocumented feature)
		% 0 - normal mode  (V12L,V34L)
		% 1 - tape mode 1  (V12M,V34M)
		% 2 - tape mode 2  (V12M,V34M)
		% 3 - tape mode 3  (V1M,V2M,V3M,V4M)
		clear tm mTMode1 mTMode2 mTMode3 mTMode4
		if exist('./mTMode.mat','file'), load mTMode, end
		if exist(irf_ssub('mTMode?',cl_id),'var')
			c_eval('tm=mTMode?;',cl_id)
		end
		if ~exist('tm','var')
			[t,data] = caa_is_get(cdb.db,start_time,dt,cl_id,'efw','FDM');
			if ~isempty(data), tm = data(5,:);
			else
				irf_log('dsrc','Cannot fetch FDM')
				cd(old_pwd)
				return
			end	       
			if tm~=tm(1)*ones(size(tm))
				irf_log('dsrc','tape mode changes during the selected time inteval')
			end
			c_eval('mTMode?=tm;',cl_id)
			if exist('./mTMode.mat','file')
				eval(irf_ssub('save -append mTMode mTMode?;',cl_id))
			else, eval(irf_ssub('save mTMode mTMode?;',cl_id))
			end
		end
		if tm<1e-30, param='10Hz'; else, param='180Hz'; end
		clear tm
	
		if start_time>toepoch([2001 07 31 00 12 29.5])&start_time<toepoch([2001 09 02 23 50 00])
			% all sc run on 180Hz filter in august 2001
			param='180Hz';
		elseif start_time>toepoch([2001 09 10 04 21 57.6])& start_time<toepoch([2001 09 15 06 30 00])
			% this needs to be investigated.... 
			param='180Hz';
		elseif start_time>toepoch([2001 07 23 00 00 00])&cl_id==2
			% 10Hz filter problem on SC2
			param='180Hz';
		end
	end
	if (start_time>toepoch([2003 9 29 0 31 0])) & (cl_id==1|cl_id==3)
		% FSW 2.4 Use P32 on SC1 and SC3
		pl={'32','34'};
		irf_log('dsrc',sprintf('            !Use p32 for sc%d',cl_id));
	elseif (start_time>toepoch([2001 12 28 03 00 00])&cl_id==1) | (start_time>toepoch([2002 07 29 09 06 59 ])&cl_id==3)
		% p1 problems on SC1 and SC3
		pl={'34'};
		irf_log('dsrc',sprintf('            !Only p34 exists for sc%d',cl_id));
	else, pl={'12','34'};
	end
	
	for i=1:length(pl)
		irf_log('dsrc',['EFW...sc' num2str(cl_id) '...Ep' pl{i} ' ' param ' filter']);
		[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, ...
		'efw', 'E', ['p' pl{i}], param, tmmode);
		if ~isempty(data)
			data = double(real(data));
			t = double(t);
			
			if do_burst
				% correct start time of the burst
				burst_f_name = irf_ssub([irf_fname(t(1),1) 'we.0?'],cl_id);
				burst_f_name = [cdb.dp '/burst/' burst_f_name];
				if exist(burst_f_name,'file')
					db = Mat_DbOpen(cdb.db);
					err_t = t(1) - c_efw_burst_chkt(db,burst_f_name);
					irf_log('dsrc',['burst start time was corrected by ' num2str(err_t) ' sec'])
					Mat_DbClose(db);
					t = t - err_t;
				else
					irf_log('dsrc','burst start time was not corrected')
				end
			end
					
			c_eval([var_name  pl{i} '=[t data];'],cl_id); clear t data;
			c_eval(['save_list=[save_list '' ' var_name pl{i} ' ''];'],cl_id);
		else
			irf_log('dsrc',irf_ssub(['No data for ' var_name pl{i}],cl_id))
		end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p - SC potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'p') | strcmp(quantity,'pburst')
	
	if strcmp(quantity,'pburst'), do_burst = 1; else do_burst = 0; end
	if do_burst 
		save_file = './mEFWburst.mat';
		tmmode='burst';
		param={'4kHz','32kHz'};
		var_name = 'wbE?p';
	else
		save_file = './mPR.mat';
		param={'10Hz'}; tmmode='lx';
	end
	
	probe_list = 1:4;
	
	% Check for p1 problems on SC1 and SC3
	if (start_time>toepoch([2001 12 28 03 00 00])&cl_id==1) | (start_time>toepoch([2002 07 29 09 06 59 ])&cl_id==3)
		probe_list = 2:4;
		p1 = [];
		irf_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
	elseif start_time>toepoch([2001 07 23 00 00 00])&cl_id==2 & ~do_burst
		probe_list = [1 2 4];
		p3 = [];
		irf_log('dsrc',sprintf('10Hz filter problem on sc%d',cl_id));
	end
	
	for j=1:length(param), for probe=probe_list;
    	irf_log('dsrc',['EFW...sc' num2str(cl_id) '...probe' num2str(probe) '->P' param{j} num2str(cl_id) 'p' num2str(probe)]);
		[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, ...
		'efw', 'E', ['p' num2str(probe)],param{j}, tmmode);
		if ~isempty(data)
			data = double(real(data));
			t = double(t);
			
			if do_burst
				% correct start time of the burst
				burst_f_name = irf_ssub([irf_fname(t(1),1) 'we.0?'],cl_id);
				burst_f_name = [cdb.dp '/burst/' burst_f_name];
				if exist(burst_f_name,'file')
					db = Mat_DbOpen(cdb.db);
					err_t = t(1) - c_efw_burst_chkt(db,burst_f_name);
					irf_log('dsrc',['burst start time was corrected by ' num2str(err_t) ' sec'])
					Mat_DbClose(db);
					t = t - err_t;
				else
					irf_log('dsrc','burst start time was not corrected')
				end
			end
			eval(irf_ssub(['p!=[t data];save_list=[save_list ''P' param{j} '?p! ''];P' param{j} '?p!=p!;'],cl_id,probe)); clear t data
		else
			eval(['p' num2str(probe) '=[];'])
		end
    end, end
	
    clear p
	if do_burst
		% make electric field
		for j=1:length(param)
			for probe=[1 3]
				if exist(irf_ssub(['P' param{j} '?p!'],cl_id,probe),'var') & ...
				exist(irf_ssub(['P' param{j} '?p!'],cl_id,probe+1),'var')
					eval(irf_ssub(['E(:,1)=P' param{j} '?p$(:,1);E(:,2)=(P' param{j} '?p$(:,2)-P' param{j} '?p!(:,2))/.088;'],cl_id,probe,probe+1));
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
	clear p1 p2 p3 p4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Phase, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'a')
	save_file = './mA.mat';
	
	% We ask for 5 sec more from each side to avoid problemos with interpolation.
	[t,data] = caa_is_get(cdb.db, start_time-5, dt+10, cl_id, 'ephemeris', 'phase');
	if ~isempty(data)
		c_eval('A?=[double(t) double(data)];',cl_id); clear t data;
		c_eval('save_list=[save_list '' A? ''];',cl_id);
	else
		irf_log('dsrc',irf_ssub('No data for A?',cl_id))
	end
	[t,data] = caa_is_get(cdb.db, start_time-5, dt+10, cl_id, 'ephemeris', 'phase_2');
	if ~isempty(data)
		c_eval('Atwo?=[double(t) double(data)];',cl_id); clear t data;
		c_eval('save_list=[save_list '' Atwo? ''];',cl_id);
	else
		irf_log('dsrc',irf_ssub('No data for Atwo?',cl_id))
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'r')
	save_file = './mR.mat';
	[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, 'ephemeris', 'position');
	if ~isempty(data)
		c_eval('R?=[double(t) double(data'')];',cl_id); clear t data;
		c_eval('save_list=[save_list '' R? ''];',cl_id);
	else
		irf_log('dsrc',irf_ssub('No data for R?',cl_id))
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'v')
	save_file = './mR.mat';
	[t,data] = caa_is_get(cdb.db, start_time, dt, cl_id, 'ephemeris', 'velocity');
	if ~isempty(data)
		c_eval('V?=[double(t) double(data'')];',cl_id); clear t data;
		c_eval('save_list=[save_list '' V? ''];',cl_id);
		% transform vector data to DSI
		if c_load('SAX?',cl_id)
			c_eval('diV?=c_gse2dsi(V?,SAX?);save_list=[save_list '' diV?''];',cl_id);
		else
			irf_log('load','must fetch spin axis orientation (option ''sax'')')
		end
	else
		irf_log('dsrc',irf_ssub('No data for V?',cl_id))
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
	disp('CONTACT STEPHAN BUCHERT!!!!!!!!!!!!!!!!!!!');

	save_file = './mB.mat';
	
	dat = c_get_bfgm(start_time + [0 dt],cl_id);
	if ~isempty(dat)
		c_eval('B?=dat;save_list=[save_list '' B?''];',cl_id);

		% transform vector data to DSI
		if size(dat,2)>2 
			if exist('./mEPH.mat','file'), eval(irf_ssub('load mEPH SAX?',cl_id)), end
			if ~exist(irf_ssub('SAX?',cl_id),'var')
				irf_log('load','must fetch spin axis orientation (option ''sax'')')
			else
				c_eval('diB?=c_gse2dsi(B?,SAX?);save_list=[save_list '' diB?''];',cl_id);
			end
		end
		clear dat
	else
		irf_log('dsrc','No data')
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CSDS PP [GSE+DSI] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'b')|strcmp(quantity,'edi')|strcmp(quantity,'vcis')|strcmp(quantity,'ncis')

	inertial_f = 0;
	r.qua = {quantity};
	switch(quantity)
	case 'b'
		r.ins = 'BPP';
		r.var = {r.ins};
	case 'edi'
		r.ins = 'EDI';
		r.var = {r.ins};
		inertial_f = 1; % EDI pp is in the inertial frame
	case 'vcis'
		r.ins = 'CIS';
		r.qua = {'vcis_p', 'vcis_h'}; % CODIF and HIA 
		r.var = {'VCp', 'VCh'};
	case 'ncis'
		r.ins = 'CIS';
		r.qua = {'ncis_p', 'ncis_h'}; % CODIF and HIA 
		r.var = {'NCp', 'NCh'};
	otherwise
		error('Check variable list')
	end

	save_file = ['./m' r.ins '.mat'];

	for i=1:length(r.qua)	
		% first try ISDAT (fast) then files
		dat = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,r.qua{i});
		if ~isempty(dat)
			if inertial_f
				c_eval(['i' r.var{i} '?=dat;save_list=[save_list '' i' r.var{i} '?''];'],cl_id);
			else
				c_eval([r.var{i} '?=dat;save_list=[save_list '' ' r.var{i} '?''];'],cl_id);
			end

			% transform vector data to DSI
			if size(dat,2)>2 
				if c_load('SAX?',cl_id)
					if inertial_f
						c_eval(['idi' r.var{i} '?=c_gse2dsi(dat,SAX?);save_list=[save_list '' idi' r.var{i} '?''];'],cl_id);
					else
						c_eval(['di' r.var{i} '?=c_gse2dsi(dat,SAX?);save_list=[save_list '' di' r.var{i} '?''];'],cl_id);
					end
				else
					irf_log('load','must fetch spin axis orientation (option ''sax'')')
				end
			end
			clear dat
		else
			irf_log('dsrc','No data')
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sax - spin axis orientation [GSE] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'sax')
	save_file='./mEPH.mat';

	% first try ISDAT (fast) then files
	lat = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slat');
	long = c_csds_read([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slong');
	if ~isempty(lat) & ~isempty(long)
		% take first point only. This is OK according to AV
		[xspin,yspin,zspin] = sph2cart(long(1,2)*pi/180,lat(1,2)*pi/180,1);
		sax = [xspin yspin zspin];
		eval(irf_ssub('SAX?=sax;save_list=[save_list '' SAX?''];',cl_id));
		clear sax lat long xspin yspin zspin
	end

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
if flag_save==1 & length(save_file)>0 & ~isempty(save_list)
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

