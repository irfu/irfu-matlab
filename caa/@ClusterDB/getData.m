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
%	p   : P{cl_id},NVps{cl_id}, P10Hz{cl_id}p{1:4} -> mP	
%			// probe potential (LX)
%
%	//// EFW internal burst////
%	eburst: wbE{cl_id}p12,34 -> mEFWburst
%			// electric fields 8kHz
%	pburst: P{4kHz,32kHz}{cl_id}p{1:4}, wbE{cl_id}p12,34 -> mEFWburst	
%			// probe potentials (4kHz,32kHz), and electric fields
%
%	//// Ephemeris ////
%	sax : SAX{cl_id} ->mEPH
%			// spin axis vector [GSE] 
%	a   : A{cl_id} -> mA	// SC phase
%	r   : R{cl_id} -> mR	// SC position
%	v   : V{cl_id} -> mR	// SC velocity
%
%	//// Other instruments ////
%	b   : BPP{cl_id},diBPP{cl_id}	->mBPP	// B FGM PP [GSE+DSI] 
%	bfgm: B{cl_id},diB{cl_id}	->mB	// B FGM** [GSE+DSI]
%		** contact Stephan Buchert
%	edi : EDI{cl_id},diEDI{cl_id}	->mEDI	// E EDI PP [GSE+DSI] 
%	ncis: NC(p,h){cl_id}			->mCIS	// N CIS PP 
%	vcis: VC(p,h){cl_id},diVC(p,h){cl_id}  ->mCIS	// V CIS PP [GSE+DSI] 
%	vce : VCE(p,h){cl_id},diVCE(p,h){cl_id} ->mCIS	// E CIS PP [GSE+DSI] 
%	wbdwf: wfWBD{cl_id} -> mWBD	// WBD waveforms E/B 
%	whip: WHIP{cl_id} -> mFDM	// Whisper pulses present +1 precceding sec 
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

% Copyright 2004 Yuri Khotyaintsev
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
		c_log('fcal',['Option ''' varargin{i} '''not recognized'])
	end
end

start_date_str = strrep(datestr(fromepoch(start_time),29),'-','');

save_file = '';
save_list = '';

old_pwd = pwd;

%Create the storage directory if it does not exist
if ~exist(cdb.sp, 'dir')
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir(cdb.sp);
	if SUCCESS, c_log('save',['Created storage directory ' cdb.sp])
	else, error(MESSAGE)
	end
end

cd(cdb.sp) %enter the storage directory
c_log('save',['Storage directory is ' cdb.sp])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e - Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'e') | strcmp(quantity,'eburst')
	
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
		if exist(av_ssub('mTMode?',cl_id),'var')
			c_eval('tm=mTMode?;',cl_id)
		end
		if ~exist('tm','var')
			[t,data] = ISGet(cdb.db,start_time,dt,cl_id,'efw','FDM');
			if ~isempty(data), tm = data(5,:);
			else, error('Cannot fetch FDM')
			end	       
			if tm~=tm(1)*ones(size(tm))
				warning('tape mode changes during the selected time inteval')
			end
			c_eval('mTMode?=tm;',cl_id)
			if exist('./mTMode.mat','file')
				eval(av_ssub('save -append mTMode mTMode?;',cl_id))
			else, eval(av_ssub('save mTMode mTMode?;',cl_id))
			end
		end
		if tm<1e-30, param='10Hz'; else, param='180Hz'; end
		clear tm
	
		if start_time>toepoch([2001 07 31 00 00 00])&start_time<toepoch([2001 09 01 00 00 00])
			% all sc run on 180Hz filter in august 2001
			param='180Hz';
		elseif start_time>toepoch([2001 09 15 04 30 00])& start_time<toepoch([2001 09 15 06 30 00])
			% this needs to be investigated.... 
			param='180Hz';
		elseif start_time>toepoch([2001 07 31 00 00 00])&cl_id==2
			% 10Hz filter problem on SC2
			param='180Hz';
		end
	end
	if (start_time>toepoch([2001 12 28 03 00 00])&cl_id==1) | (start_time>toepoch([2002 07 29 09 06 59 ])&cl_id==3)
		% p1 problems on SC1 and SC3
		pl={'34'};
		c_log('dsrc',sprintf('            !Only p34 exists for sc%d',cl_id));
	else, pl={'12','34'};
	end
	for i=1:length(pl)
		c_log('dsrc',['EFW...sc' num2str(cl_id) '...Ep' pl{i} ' ' param ' filter']);
		[t,data] = ISGet(cdb.db, start_time, dt, cl_id, ...
		'efw', 'E', ['p' pl{i}], param, tmmode);
		if ~isempty(data)
			data = double(real(data));
			t = double(t);
			
			if do_burst
				% correct start time of the burst
				burst_f_name = av_ssub([makeFName(t(1),1) 'we.0?'],cl_id);
				burst_f_name = [cdb.dp '/burst/' burst_f_name];
				if exist(burst_f_name,'file')
					db = Mat_DbOpen(cdb.db);
					err_t = t(1) - checkbursttime(db,burst_f_name);
					c_log('dsrc',['burst start time was corrected by ' num2str(err_t) ' sec'])
					Mat_DbClose(db);
					t = t - err_t;
				else
					c_log('dsrc','burst start time was not corrected')
				end
			end
					
			c_eval([var_name  pl{i} '=[t data];'],cl_id); clear t data;
			c_eval(['save_list=[save_list '' ' var_name pl{i} ' ''];'],cl_id);
		else
			c_log('dsrc',av_ssub(['No data for ' var_name pl{i}],cl_id))
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
		save_file = './mP.mat';
		param={'10Hz'}; tmmode='lx';
	end
	
	probe_list = 1:4;
	
	% Check for p1 problems on SC1 and SC3
	if (start_time>toepoch([2001 12 28 03 00 00])&cl_id==1) | (start_time>toepoch([2002 07 29 09 06 59 ])&cl_id==3)
		probe_list = 2:4;
		p1 = [];
		c_log('dsrc',sprintf('p1 is BAD on sc%d',cl_id));
	end
	
	for j=1:length(param), for probe=probe_list;
    	c_log('dsrc',['EFW...sc' num2str(cl_id) '...probe' num2str(probe) '->P' param{j} num2str(cl_id) 'p' num2str(probe)]);
		[t,data] = ISGet(cdb.db, start_time, dt, cl_id, ...
		'efw', 'E', ['p' num2str(probe)],param{j}, tmmode);
		if ~isempty(data)
			data = double(real(data));
			t = double(t);
			
			if do_burst
				% correct start time of the burst
				burst_f_name = av_ssub([makeFName(t(1),1) 'we.0?'],cl_id);
				burst_f_name = [cdb.dp '/burst/' burst_f_name];
				if exist(burst_f_name,'file')
					db = Mat_DbOpen(cdb.db);
					err_t = t(1) - checkbursttime(db,burst_f_name);
					c_log('dsrc',['burst start time was corrected by ' num2str(err_t) ' sec'])
					Mat_DbClose(db);
					t = t - err_t;
				else
					c_log('dsrc','burst start time was not corrected')
				end
			end
			eval(av_ssub(['p!=[t data];save_list=[save_list '' P' param{j} '?p!''];P' param{j} '?p!=p!;'],cl_id,probe)); clear t data
		else
			eval(['p' num2str(probe) '=[];'])
		end
    end, end
	
    clear p
	if do_burst
		% make electric field
		for j=1:length(param)
			for probe=[1 3]
				if exist(av_ssub(['P' param{j} '?p!'],cl_id,probe),'var') & ...
				exist(av_ssub(['P' param{j} '?p!'],cl_id,probe+1),'var')
					eval(av_ssub(['E(:,1)=P' param{j} '?p$(:,1);E(:,2)=(P' param{j} '?p$(:,2)-P' param{j} '?p!(:,2))*12.5;'],cl_id,probe,probe+1));
					vn = [var_name num2str(probe) num2str(probe+1)];
					if exist(av_ssub(vn,cl_id),'var')
						c_eval(['tmpE=' vn ';'],cl_id)
						if tmpE(1,1) > E(1,1)
							tmpE(:,end+1:end+size(E,1)) = E;
							E = tmpE;
						else
							E(:,end+1:end+size(tmpE,1)) = tmpE;
						end
						clear tmpE
					end
					c_eval([vn '=E;save_list=[save_list '' ' vn  ' ''];'],cl_id);
					clear E
				end
			end
		end
	else
		if size(p1)==size(p2)&size(p1)==size(p3)&size(p1)==size(p4)&size(p1)~=[0 0]&cl_id~=2,  % sc2 has often problems with p3
		p=[p1(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4];
		elseif size(p1)==size(p2)&size(p1)~=[0 0]
				p=[p1(:,1) (p1(:,2)+p2(:,2))/2];
		elseif size(p3)==size(p4)&cl_id~=2
			p=[p3(:,1) (p3(:,2)+p4(:,2))/2];
		else,
			p=p4;
		end
		
		c_eval(['P' param{1} '?=p;save_list=[save_list '' P' param{1} '? ''];'],cl_id);
		c_eval('P?=p;NVps?=c_n_Vps(p);NVps?(:,end+1)=p(:,2); save_list=[save_list '' P? NVps?''];',cl_id)
		clear p
	end
	clear p1 p2 p3 p4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Phase, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'a')
	save_file = './mA.mat';
	[t,data] = ISGet(cdb.db, start_time, dt, cl_id, 'ephemeris', 'phase');
	if ~isempty(data)
		c_eval('A?=[double(t) double(data)];',cl_id); clear t data;
		c_eval('save_list=[save_list '' A? ''];',cl_id);
	else
		c_log('dsrc',av_ssub('No data for A?',cl_id))
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'r')
	save_file = './mR.mat';
	[t,data] = ISGet(cdb.db, start_time, dt, cl_id, 'ephemeris', 'position');
	if ~isempty(data)
		c_eval('R?=[double(t) double(data'')];',cl_id); clear t data;
		c_eval('save_list=[save_list '' R? ''];',cl_id);
	else
		c_log('dsrc',av_ssub('No data for R?',cl_id))
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'v')
	save_file = './mR.mat';
	[t,data] = ISGet(cdb.db, start_time, dt, cl_id, 'ephemeris', 'velocity');
	if ~isempty(data)
		c_eval('V?=[double(t) double(data'')];',cl_id); clear t data;
		c_eval('save_list=[save_list '' V? ''];',cl_id);
	else
		c_log('dsrc',av_ssub('No data for V?',cl_id))
	end

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
			if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
			if ~exist(av_ssub('SAX?',cl_id),'var')
				c_log('load','must fetch spin axis orientation (option ''sax'')')
			else
				c_eval('diB?=c_gse2dsi(B?,SAX?);save_list=[save_list '' diB?''];',cl_id);
			end
		end
		clear dat
	else
		c_log('dsrc','No data')
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CSDS PP [GSE+DSI] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'b')|strcmp(quantity,'edi')|strcmp(quantity,'vcis')|strcmp(quantity,'ncis')

	r.qua = {quantity};
	switch(quantity)
	case 'b'
		r.ins = 'BPP';
		r.var = {r.ins};
	case 'edi'
		r.ins = 'EDI';
		r.var = {r.ins};
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
		dat = readCSDS([cdb.db '|' cdb.dp],start_time,dt,cl_id,r.qua{i});
		if ~isempty(dat)
			eval(av_ssub([r.var{i} '?=dat;save_list=[save_list '' ' r.var{i} '?''];'],cl_id));

			% transform vector data to DSI
			if size(dat,2)>2 
				if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
				if ~exist(av_ssub('SAX?',cl_id),'var')
					c_log('load','must fetch spin axis orientation (option ''sax'')')
				else
					eval(av_ssub(['di' r.var{i} '?=c_gse2dsc(' r.var{i} '?,SAX?);di' r.var{i} '?(:,3:4)=-di' r.var{i} '?(:,3:4);save_list=[save_list '' di' r.var{i} '?''];'],cl_id));
				end
			end
			clear dat
		else
			c_log('dsrc','No data')
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sax - spin axis orientation [GSE] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'sax')
	save_file='./mEPH.mat';

	% first try ISDAT (fast) then files
	lat = readCSDS([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slat');
	long = readCSDS([cdb.db '|' cdb.dp],start_time,dt,cl_id,'slong');
	if ~isempty(lat) & ~isempty(long)
		% take first point only. This is OK according to AV
		[xspin,yspin,zspin] = sph2cart(long(1,2)*pi/180,lat(1,2)*pi/180,1);
		sax = [xspin yspin zspin];
		eval(av_ssub('SAX?=sax;save_list=[save_list '' SAX?''];',cl_id));
		clear sax lat long xspin yspin zspin
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vce - E CIS PP [GSE+DSI] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'vce')
	save_file='./mCIS.mat';

	if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
	if ~exist('./mCIS.mat','file')
		c_log('load','Please run ''vcis'' first (mCIS missing)')
		data = [];
		return
	end
	if ~exist('./mBPP.mat','file')
		c_log('load','Please run ''n'' first (mBPP missing)')
		data = [];
		return
	end

	CIS=load('mCIS');
	eval(av_ssub('load mBPP BPP?;b=BPP?; clear BPP?',cl_id))

	var = {'VCp', 'VCh'};
	varo = {'VCEp', 'VCEh'};
	for i=1:length(var)
	
		eval(av_ssub(['if isfield(CIS,''' var{i} '?''); v=CIS.' var{i} '?; else, v=[]; end; clear ' var{i}], cl_id));
		if ~isempty(v)
			evxb=av_t_appl(av_cross(v,b),'*(-1e-3)');
			eval(av_ssub([varo{i} '?=evxb;save_list=[save_list '' ' varo{i} '?'']; clear evxb'],cl_id));
			if ~exist(av_ssub('SAX?',cl_id),'var')
				c_log('load','must fetch spin axis orientation (option ''sax'')')
			else
				eval(av_ssub(['di' varo{i} '?=c_gse2dsc(' varo{i} '?,SAX?);di' varo{i} '?(:,3:4)=-di' varo{i} '?(:,3:4);save_list=[save_list '' di' varo{i} '?''];'],cl_id));
			end
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whip - Whisper present.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'whip')
	save_file = './mFDM.mat';
	[t_s,t_e,fdm_r]=createEFWModeTableFDM(cdb.db, start_time, dt, cl_id, 'r');
	ii = find(fdm_r==1);
	
	if ~isempty(ii)
		% add 1 sec from both sides
		t_s = t_s(ii) - 1;
		t_e = t_e(ii);
		c_eval('WHIP?=[double(t_s)'' double(t_e)''];',cl_id); 
		c_eval('save_list=[save_list '' WHIP? ''];',cl_id);
	else
		c_log('dsrc','No data')
	end
	clear t_s t_e fdm_r ii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wbdwf - WBD waveforms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'wbdwf')
	save_file = './mWBD.mat';
	try
		wf = readWBD(start_time, dt, cl_id);
		c_eval('wfWBD?=wf;',cl_id); 
		c_eval('save_list=[save_list '' wfWBD? ''];',cl_id);
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else, error('caa:noSuchQuantity','Quantity ''%s'' unknown',quantity)
end %main QUANTITY

% saving
% If flag_save is set, save variables to specified file
if flag_save==1 & length(save_file)>0 & ~isempty(save_list)
	c_log('save',[save_list ' -> ' save_file])
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

