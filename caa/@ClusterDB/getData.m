function data = getData(cdb,start_time,dt,cl_id,quantity,varargin)
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
%	e : wE{cl_id}p12, wE{cl_id}p34 -> mER // electric fields
%	a : A{cl_id} -> mA // phase
%
%	options - one of the following:
%	not yet implemented
%
% $Revision$  $Date$
%
% see also C_GET, TOEPOCH

% Copyright 2004 Yuri Khotyaintsev
% Parts of the code are (c) Andris Vaivads

error(nargchk(5,15,nargin))
if nargin > 5, property_argin = varargin; end

flag_save = 1; % change this!

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
% aux data - Phase, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(quantity,'a')
	save_file = './mA.mat';
	[t,data] = ISGet(cdb.db, start_time, dt, cl_id, 'ephemeris', 'phase');
	eval(av_ssub('A?=[double(t) double(data)];',cl_id)); clear t data;
	eval(av_ssub('save_list=[save_list '' A? ''];',cl_id));
end %main QUANTITY

% saving
% If flag_save is set, save variables to specified file
if flag_save==1 & length(save_file)>0 & ~isempty(save_list)
	if exist(save_file,'file')
		eval(['save -append ' save_file ' ' save_list]);
	else
		eval(['save ' save_file ' ' save_list]);
	end
end

% prepare the output
if nargout > 0 & ~isempty(save_list)
	sl = tokenize(save_list);
	data = {sl};
	for i=1:length(sl)
		eval(['data{i+1}={' char(sl{i}) '};'])
	end
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

% reset errors, otherwise try/catch fails
lasterr('')
try
	dbase = Mat_DbOpen(db_s);
	disp(lasterr)

	[t, d] = isGetDataLite(dbase, st, dt, ...
	'Cluster',num2str(cli),ins,sig,sen,cha,par);

	if exist('dbase'), Mat_DbClose(dbase); end

catch
	disp(sprintf('Error getting data from ISDAT server %s: %s',...
	db_s,lasterr))
	t = []; d = [];
end
