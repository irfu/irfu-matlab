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
%	e   : wE{cl_id}p12,34 -> mER
%			// electric fields (HX)
%	p   : P{cl_id},NVps{cl_id}, P10Hz{cl_id}p{1:4} -> mP	
%			// sc potential (LX)
%	sax : SAX{cl_id} ->mEPH
%			// spin axis vector [GSE] 
%	a   : A{cl_id} -> mA	// phase
%	b   : BPP{cl_id},diBPP{cl_id} ->mBPP	// B FGM PP [GSE+DSI] 
%	edi : EDI{cl_id},diEDI{cl_id} ->mEDI	// E EDI PP [GSE+DSI] 
%	vcis: VC(p,h){cl_id},diVC(p,h){cl_id}  ->mCIS	// V CIS PP [GSE+DSI] 
%	vce: VCE(p,h){cl_id},diVCE(p,h){cl_id} ->mCIS	// E CIS PP [GSE+DSI] 
%
%	options - one of the following:
%	nosave : do no save on disk
%
% $Revision$  $Date$
%
% see also C_GET, TOEPOCH

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
		disp(['Option ''' varargin{i} '''not recognized'])
	end
end

start_date_str = strrep(datestr(fromepoch(start_time),29),'-','');

save_file = '';
save_list = '';

old_pwd = pwd;
cd(cdb.sp) %enter the storage directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e - Electric field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(quantity,'e')
	save_file = './mER.mat';
	tmmode='hx';
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
		eval(av_ssub('tm=mTMode?;',cl_id))
	end
	if ~exist('tm','var')
		[t,data] = ISGet(cdb.db,start_time,dt,cl_id,'efw','FDM');
		if ~isempty(data), tm = data(5,:);
		else, error('Cannot fetch FDM')
		end	       
		if tm~=tm(1)*ones(size(tm))
			warning('tape mode changes during the selected time inteval')
		end
		eval(av_ssub('mTMode?=tm;',cl_id))
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
	elseif start_time>toepoch([2001 07 31 00 00 00])&cl_id==2
		% 10Hz filter problem on SC2
		param='180Hz';
	end
	if (start_time>toepoch([2001 12 28 03 00 00])&cl_id==1) | (start_time>toepoch([2002 07 29 09 06 59 ])&cl_id==3)
		% p1 problems on SC1 and SC3
		pl={'34'};
		disp(sprintf('            !Only p34 exists for sc%d',cl_id));
	else, pl={'12','34'};
	end
	for i=1:length(pl)
		disp(['EFW...sc' num2str(cl_id) '...Ep' pl{i} ' ' param ' filter']);
		[t,data] = ISGet(cdb.db, start_time, dt, cl_id, ...
		'efw', 'E', ['p' pl{i}], param, tmmode);
		data = double(real(data));
		t = double(t);
		eval(av_ssub(['wE?p'  pl{i} '=[t data];'],cl_id)); clear t data;
		eval(av_ssub(['save_list=[save_list '' wE?p' pl{i} ' ''];'],cl_id));
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p - SC potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'p')
	save_file = './mP.mat';
	param='10Hz'; tmmode='lx';
	for probe=1:4;
    	disp(['EFW...sc' num2str(cl_id) '...probe' num2str(probe) '->P' param num2str(cl_id) 'p' num2str(probe)]);
		[t,data] = ISGet(cdb.db, start_time, dt, cl_id, ...
		'efw', 'E', ['p' num2str(probe)],param, tmmode);
		if ~isempty(data)
			eval(av_ssub(['p!=[double(t) double(real(data))];save_list=[save_list '' P' param '?p!''];P' param '?p!=p!;'],cl_id,probe)); clear t data
		else
			eval(['p' num2str(probe) '=[];'])
		end
    end
    clear p
	if size(p1)==size(p2)&size(p1)==size(p3)&size(p1)==size(p4)&size(p1)~=[0 0]&cl_id~=2,  % sc2 has often problems with p3
    	p=[p1(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4];
    elseif size(p1)==size(p2)&size(p1)~=[0 0]
    	p=[p1(:,1) (p1(:,2)+p2(:,2))/2];
    elseif size(p3)==size(p4)&cl_id~=2
    	p=[p3(:,1) (p3(:,2)+p4(:,2))/2];
    else,
    	p=p4;
    end
    eval(av_ssub(['P' param '?=p;save_list=[save_list '' P' param '? ''];'],cl_id));
	eval(av_ssub('P?=p;NVps?=c_n_Vps(p);NVps?(:,end+1)=p(:,2); save_list=[save_list '' P? NVps?''];',cl_id))
	clear p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux data - Phase, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'a')
	save_file = './mA.mat';
	[t,data] = ISGet(cdb.db, start_time, dt, cl_id, 'ephemeris', 'phase');
	if ~isempty(data)
		eval(av_ssub('A?=[double(t) double(data)];',cl_id)); clear t data;
		eval(av_ssub('save_list=[save_list '' A? ''];',cl_id));
	else
		warning('caa:noData',av_ssub('No data for A?',cl_id))
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CSDS PP [GSE+DSI] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'b')|strcmp(quantity,'edi')|strcmp(quantity,'vcis')
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
	otherwise
		error('Check variable list')
	end

	save_file = ['./m' r.ins '.mat'];

	for i=1:length(r.qua)	
		% first try ISDAT (fast) then files
		dat = readCSDS([cdb.db '|' cdb.dp],start_time,dt,cl_id,r.qua{i});
		if ~isempty(dat)
			eval(av_ssub([r.var{i} '?=dat;save_list=[save_list '' ' r.var{i} '?''];'],cl_id));
			if exist('./mEPH.mat','file'), eval(av_ssub('load mEPH SAX?',cl_id)), end
			if ~exist(av_ssub('SAX?',cl_id),'var')
				warning('must fetch spin axis orientation (option ''sax'')')
			else
				eval(av_ssub(['di' r.var{i} '?=c_gse2dsc(' r.var{i} '?,SAX?);di' r.var{i} '?(:,3:4)=-di' r.var{i} '?(:,3:4);save_list=[save_list '' di' r.var{i} '?''];'],cl_id));
			end
			clear dat
		else
			warning('caa:noData','No data')
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
		warning('Please run ''vcis'' first (mCIS missing)')
		data = [];
		return
	end
	if ~exist('./mBPP.mat','file')
		warning('Please run ''n'' first (mBPP missing)')
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
				warning('must fetch spin axis orientation (option ''sax'')')
			else
				eval(av_ssub(['di' varo{i} '?=c_gse2dsc(' varo{i} '?,SAX?);di' varo{i} '?(:,3:4)=-di' varo{i} '?(:,3:4);save_list=[save_list '' di' varo{i} '?''];'],cl_id));
			end
		end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else, error('caa:noSuchQuantity','Quantity ''%s'' unknown',quantity)
end %main QUANTITY

% saving
% If flag_save is set, save variables to specified file
if flag_save==1 & length(save_file)>0 & ~isempty(save_list)
	disp([save_list ' -> ' save_file])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help functions - ISGet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,d] = ISGet(db_s,st,dt,cli,ins,sig,sen,cha,par)
% Help wrapper around isGetDataLite
% Last 3 input arguments can be skipped

if nargin < 9, par = ' '; end
if nargin < 8, cha = ' '; end
if nargin < 7, sen = ' '; end

p = tokenize(db_s,'|');

for i=1:length(p)
	% reset errors, otherwise try/catch fails
	lasterr('')
	try
		dbase = Mat_DbOpen(p{i});
		disp(lasterr)

		[t, d] = isGetDataLite(dbase, st, dt, ...
		'Cluster',num2str(cli),ins,sig,sen,cha,par);

		if exist('dbase'), Mat_DbClose(dbase), clear dbase, end

		if ~isempty(d), return, end

	catch
		warning('ISDAT:getData',...
		'Error getting data from ISDAT server %s: %s', p{i},lasterr)
		t = []; d = [];
	end
end
