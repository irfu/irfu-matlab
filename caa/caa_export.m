function status = caa_export(lev,caa_vs,cl_id,QUALITY,DATA_VERSION,sp,st,dt)
%CAA_EXPORT  export data to CAA CEF files
%
% STATUS = caa_export(data_level,caa_vs,cl_id,[QUALITY,DATA_VERSION,sp,st,dt])
%
% CEF file is written to a current directory. Data is supposed to be also 
% there, if SP is not specified.
%
% QUALITY - number 0..4 (3 is good for publication, subject to PI approval)
% DATA_VERSION - string, for example '01'
%
% STATUS = 0 means everything went OK
%
% See also C_EXPORT_ASCII
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

status = 0;

% This must be changed when we do any major changes to our processing software
EFW_DATASET_VERSION = '2';

if nargin<8, st = []; dt=[]; end
if nargin<6, sp='.'; end
if nargin<5, DATA_VERSION = '02'; end
if nargin<4, QUALITY = 3; end % Good for publication, subject to PI approval
if cl_id<=0 || cl_id>4, error('CL_ID must be 1..4'), end
if lev<1 || lev>3, error('LEV must be 1,2 or 3'), end

DATASET_DESCRIPTION_PREFIX = '';
EOR_MARKER = '$';
FILL_VAL = -1.0E9;
PROCESSING_LEVEL='Calibrated';

old_pwd = pwd;
cd(sp)

if lev==1
	if regexp(caa_vs,'^P(1|2|3|4|12|32|34)?$')
		id = str2double(caa_vs(2:end));
		if id <=4, vs = irf_ssub('P10Hz?p!',cl_id,id);
        else vs = irf_ssub('wE?p!',cl_id,id);
		end
		v_size = 1;
    else error('Must be P(1|2|3|4|12|32|34)')
	end
else
	switch caa_vs
	case 'P'
		if lev==2
			vs = irf_ssub('P?',cl_id);
			v_size = 1;
		else
			vs = irf_ssub('Ps?',cl_id);
			v_size = 1;
		end
		DATASET_DESCRIPTION_PREFIX = 'negative of the ';
	case 'E'
		if lev==2
			vs = irf_ssub('diE?p1234',cl_id);
			v_size = 1;
			sfit_probe = caa_sfit_probe(cl_id);
		elseif lev==3
			sfit_probe = caa_sfit_probe(cl_id);
			vs = irf_ssub('diEs?p!',cl_id,sfit_probe);
			irf_log('proc',sprintf('using p%d',sfit_probe))
			v_size = 2;
		else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'DER'
		sfit_probe = caa_sfit_probe(cl_id);
		vs = irf_ssub('Dadc?p!',cl_id,sfit_probe);
		irf_log('proc',sprintf('using p%d',sfit_probe))
		v_size = 1;
	otherwise
		error('unknown variable')
	end
end

% Load data
[ok,data] = c_load(vs);
if ~ok || isempty(data)
	irf_log('load', ['No ' vs])
	cd(old_pwd)
	return
end
d_info = []; ok = 0;
try
	[ok, d_info] = c_load([vs '_info'],'var');
catch
	d_info = []; ok = 0;
end

if ~ok || isempty(d_info), dsc = c_desc(vs);
else dsc = c_desc(vs,d_info);
end

if lev==3
	TIME_RESOLUTION = 4;
	MIN_TIME_RESOLUTION = TIME_RESOLUTION;
	MAX_TIME_RESOLUTION = TIME_RESOLUTION;
elseif (lev==1 && ~isempty(regexp(caa_vs,'^P(1|2|3|4)?$','once'))) || ...
		(lev==2 && strcmp(caa_vs,'P'))
	TIME_RESOLUTION = 1/5;
	MIN_TIME_RESOLUTION = TIME_RESOLUTION;
	MAX_TIME_RESOLUTION = TIME_RESOLUTION;
elseif (lev==1 && ~isempty(regexp(caa_vs,'^P(12|32|34)?$','once'))) || ...
		(lev==2 && (strcmp(caa_vs,'E') || strcmp(caa_vs,'EF')))
	fs = c_efw_fsample(data,'hx');
	if ~fs, error('cannot determine time resolution'), end
	TIME_RESOLUTION = 1/fs;
	MIN_TIME_RESOLUTION = 1/25;
	MAX_TIME_RESOLUTION = 1/450;
end

% Make subinterval
if ~isempty(st) && ~isempty(dt)
	t_int = st + [0 dt];
	irf_log('save', sprintf('%s : %s -- %s',...
			vs, epoch2iso(t_int(1),1), epoch2iso(t_int(2),1)))
	data = irf_tlim(data,t_int);
	if isempty(data)
		irf_log('save', 'Saving empty subinterval')
	end 
else
	[iso_ts,dtint] = caa_read_interval;
	if isempty(iso_ts)
		t_int = data([1 end],1);
	else
		t_int(1) = iso2epoch(iso_ts);
		t_int(2) = t_int(1) + dtint;
	end
	clear iso_ts dtint
end

% Do magic on E-field
if strcmp(caa_vs,'E')
	
	% We export only X and Y, no need to export zeroes in Ez.
	dsc.size(1) = 2;
	
	% We check if this full res E is from coming from two probe pairs
	if lev==2 && ~(strcmp(dsc.sen,'1234') || strcmp(dsc.sen,'3234')) && QUALITY>1
		irf_log('save','This is not a full E, setting QUALITY=1!')
		QUALITY = 1;
	end
	
	% Remove wakes
	problems = 'wake'; %#ok<NASGU>
	signal = data; %#ok<NASGU>
	probe = sfit_probe; %#ok<NASGU>
	remove_problems
	data = res; %#ok<NODEF>
	clear res signal problems
	
	% Correct offsets
	if ~isempty(data)
		dsiof = c_ctl(cl_id,'dsiof');
		if isempty(dsiof)
			[ok,Ps,msg] = c_load('Ps?',cl_id);
			if ~ok, irf_log('load',msg), end
			[dsiof_def, dam_def] = c_efw_dsi_off(t_int(1),cl_id,Ps);

			[ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end
			[ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end
			
			if ok1 || ok2, irf_log('calb','Using saved DSI offsets')
			else irf_log('calb','Using default DSI offsets')
			end 
			clear dsiof_def dam_def
		else
			Ddsi = dsiof(1); Damp = dsiof(2);
			irf_log('calb','Using user specified DSI offsets')
		end
		clear dsiof
	
		data = caa_corof_dsi(data,Ddsi,Damp);
		
		if length(Ddsi) == 1
			dsc.com = sprintf('ISR2 offsets: dEx=%1.2f dEy=%1.2f dAmp=%1.2f',...
				real(Ddsi(1)),imag(Ddsi(1)),Damp);
		else
			dsc.com = 'ISR2 offsets';
			for in = 1:size(Ddsi,1)
				dsc.com = [dsc.com sprintf(' %s: dEx=%1.2f dEy=%1.2f,',...
					epoch2iso(Ddsi(in,1),1),real(Ddsi(in,2)),imag(Ddsi(in,2)))];
			end
			dsc.com = [dsc.com sprintf(' dAmp=%1.2f',Damp)];
		end
		clear Ddsi Damp
		irf_log('calb',dsc.com)
		dsc.com = [dsc.com '. Probes: ' dsc.sen];
		
		% Remove Ez, which is zero
		if lev==3, data = data(:,[1:3 5]);
        else data = data(:,1:3);
		end
	end
	
	[ok,Del] = c_load('D?p12p34',cl_id); if ~ok, Del = [0 0]; end
	if ~isreal(Del)
		Del = imag(Del);
		% offset is applied to p34
		if strcmp(dsc.sen,'12') || strcmp(dsc.sen,'32'), Del = [0 0]; end
		dsc.com = sprintf('%s. p%s offset (ISR2): dEx=%1.2f dEy=%1.2f',...
			dsc.com, '34', Del(1), Del(2));
	else
		% offset is applied to p12/32
		if strcmp(dsc.sen,'34'), Del = [0 0]; end
		dsc.com = sprintf('%s. p%s offset (ISR2): dEx=%1.2f dEy=%1.2f',...
			dsc.com, dsc.sen(1:2), Del(1), Del(2));
	end
	
	dsc.frv = {'Observatory'};
	if v_size>1, for j=2:v_size, dsc.frv = [dsc.frv {''}]; end, end
	
elseif strcmp(caa_vs,'EF')
	% We check if this full res E is from coming from two probe pairs
	if lev==2 && ~(strcmp(dsc.sen,'1234') || strcmp(dsc.sen,'3234')) && QUALITY>1
		irf_log('save','This is not a full E, setting QUALITY=1!')
		QUALITY = 1;
	end
	
	if ~isempty(data)
		% Remove Ez, which is zero
		data = data(:,1:3);
	end
	
	dsc.frv = {'Observatory'};
	
elseif lev==1 && ~isempty(regexp(caa_vs,'^P(12|32|34)?$','once'))
	if ~isempty(data)
		% convert mV/m back to V
		if id==32, data(:,2) = data(:,2)*.0622;
        else data(:,2) = data(:,2)*.088;
		end
	end
	
	dsc.units = {'V'};
	dsc.cs = {'na'};
	dsc.si_conv = {''};
	dsc.field_name = {['Potential difference measured between probes '...
		dsc.sen(1) ' and ' dsc.sen(2)]};
	dsc.ent = {'Instrument'};
	dsc.prop = {'Probe_Potential'};
	dsc.fluc = {'Waveform'};
end

cd(old_pwd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write to file
ext_s = '.cef';
% We have special names for CAA
DATASET_ID = irf_ssub(['C?_CP_EFW_L' num2str(lev) '_' caa_vs],cl_id);
file_name = ...
	[DATASET_ID '__' irf_fname(t_int,2) '_V' DATA_VERSION];

buf = '';

buf = sprintf('%s%s',buf,'!-------------------- CEF ASCII FILE --------------------|\n');
nnow = now;
buf = sprintf('%s%s',buf,['! created on ' datestr(nnow) '\n']);
buf = sprintf('%s%s',buf,'!--------------------------------------------------------|\n');
buf = sprintf('%s%s',buf,['FILE_NAME = "' file_name ext_s '"\n']);
buf = sprintf('%s%s',buf,'FILE_FORMAT_VERSION = "CEF-2.0"\n');
buf = sprintf('%s%s',buf,['END_OF_RECORD_MARKER = "' EOR_MARKER '"\n']);
buf = sprintf('%s%s',buf,'include = "CL_CH_MISSION.ceh"\n');
buf = sprintf('%s%s',buf,irf_ssub('include = "C?_CH_OBS.ceh"\n',cl_id));
buf = sprintf('%s%s',buf,'include = "CL_CH_EFW_EXP.ceh"\n');
buf = pmeta(buf,'FILE_TYPE','cef');
buf = pmeta(buf,'DATA_TYPE','CP');
buf = pmeta(buf,'INSTRUMENT_NAME','EFW?',cl_id);
buf = pmeta(buf,'INSTRUMENT_DESCRIPTION','EFW Experiment on Cluster C?',cl_id);
buf = pmeta(buf,'INSTRUMENT_CAVEATS','*C?_CQ_EFW_CAVEATS',cl_id);
buf = pmeta(buf,'DATASET_ID',DATASET_ID,cl_id);
buf = pmeta(buf,'DATASET_TITLE',dsc.field_name{1});
buf = pmeta(buf,'DATASET_DESCRIPTION',...
	{'This dataset contains measurements of the', ...
	[DATASET_DESCRIPTION_PREFIX dsc.field_name{1}],... 
	irf_ssub('from the EFW experiment on the Cluster C? spacecraft',cl_id)});
buf = pmeta(buf,'DATASET_VERSION',EFW_DATASET_VERSION);
buf = pmeta(buf,'TIME_RESOLUTION',TIME_RESOLUTION);
buf = pmeta(buf,'MIN_TIME_RESOLUTION',MIN_TIME_RESOLUTION);
buf = pmeta(buf,'MAX_TIME_RESOLUTION',MAX_TIME_RESOLUTION);
buf = pmeta(buf,'PROCESSING_LEVEL',PROCESSING_LEVEL);
buf = pmeta(buf,'DATASET_CAVEATS',['*C?_CQ_EFW_' caa_vs],cl_id);
buf = pmeta(buf,'LOGICAL_FILE_ID',file_name);
buf = pmeta(buf,'VERSION_NUMBER',DATA_VERSION);
buf = sprintf('%s%s',buf,'START_META     =   FILE_TIME_SPAN\n');
buf = sprintf('%s%s',buf,'   VALUE_TYPE  =   ISO_TIME_RANGE\n');
buf = sprintf('%s%s',buf,...
['   ENTRY       =   ' epoch2iso(t_int(1)) '/' epoch2iso(t_int(2)) '\n']);
buf = sprintf('%s%s',buf,'END_META       =   FILE_TIME_SPAN\n');
buf = sprintf('%s%s',buf,'START_META     =   GENERATION_DATE\n');
buf = sprintf('%s%s',buf,'   VALUE_TYPE  =   ISO_TIME\n');
buf = sprintf('%s%s',buf,['   ENTRY       =   ' epoch2iso(date2epoch(nnow)) '\n']);
buf = sprintf('%s%s',buf,'END_META       =   GENERATION_DATE\n');

buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,'!                   Variables                         !\n');
buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,['START_VARIABLE    = time_tags__' DATASET_ID '\n']);
buf = sprintf('%s%s',buf,'  VALUE_TYPE      = ISO_TIME\n');
buf = sprintf('%s%s',buf,['  DELTA_PLUS      = ' num2str(TIME_RESOLUTION/2) '\n']);
buf = sprintf('%s%s',buf,['  DELTA_MINUS     = ' num2str(TIME_RESOLUTION/2) '\n']);
buf = sprintf('%s%s',buf, '  FILLVAL         = 9999-12-31T23:59:59Z\n');
buf = sprintf('%s%s',buf,'  LABLAXIS        = "UT"\n');
buf = sprintf('%s%s',buf,'  FIELDNAM        = "Universal Time"\n');
buf = sprintf('%s%s',buf,['END_VARIABLE      = time_tags__' DATASET_ID '\n!\n']);

for j=1:v_size
	buf = sprintf('%s%s',buf,['START_VARIABLE      = ' dsc.name{j} '__' DATASET_ID '\n']);
	buf = sprintf('%s%s',buf,'  PARAMETER_TYPE    = "Data"\n');
	buf = sprintf('%s%s',buf,['  SIZES             = ' num2str(dsc.size(j)) '\n']);
	buf = sprintf('%s%s',buf,'  VALUE_TYPE        = FLOAT\n');
	buf = sprintf('%s%s',buf,['  ENTITY            = "' dsc.ent{j} '"\n']);
	buf = sprintf('%s%s',buf,['  PROPERTY          = "' dsc.prop{j} '"\n']);
	if ~isempty(dsc.fluc{j})
		buf = sprintf('%s%s',buf,['  FLUCTUATIONS      = "' dsc.fluc{j} '"\n']);
	end
	buf = sprintf('%s%s',buf,['  CATDESC           = "' dsc.field_name{j} '"\n']);
	buf = sprintf('%s%s',buf,['  FIELDNAM          = "' dsc.field_name{j} '"\n']);
	if ~isempty(dsc.si_conv{j})
		buf = sprintf('%s%s',buf,['  SI_CONVERSION     = "' dsc.si_conv{j} '"\n']);
	else
		buf = sprintf('%s%s',buf,['  SI_CONVERSION     = "1>' dsc.units{j} '"\n']);
	end
	buf = sprintf('%s%s',buf,['  UNITS             = "' dsc.units{j} '"\n']);
	buf = sprintf('%s%s',buf,['  FILLVAL           = ' num2str(FILL_VAL,'%8.3f') '\n']);
	buf = sprintf('%s%s',buf,['  QUALITY           = ' num2str(QUALITY) '\n']);
	buf = sprintf('%s%s',buf,'  SIGNIFICANT_DIGITS= 6 \n');
	if ~isempty(dsc.com) && j==1
		buf = sprintf('%s%s',buf,['  PARAMETER_CAVEATS = "' dsc.com '"\n']);
	end
	if ~strcmp(dsc.cs{j},'na')
		buf = sprintf('%s%s',buf,['  COORDINATE_SYSTEM = "' dsc.cs{j} '"\n']); 
	end
	if isfield(dsc,'frv')
		if ~isempty(dsc.frv{j})
			buf = sprintf('%s%s',buf,['  FRAME_VELOCITY    = "' dsc.frv{j} '"\n']);
		end
	end
	if dsc.size(j) > 1
		buf = sprintf('%s%s',buf,['  REPRESENTATION_1  = ' dsc.rep_1{j} '\n']);
		buf = sprintf('%s%s',buf,['  LABEL_1           = ' dsc.label_1{j} '\n']);
	end
	buf = sprintf('%s%s',buf,['  LABLAXIS          = "' dsc.labels{j} '"\n']);
	buf = sprintf('%s%s',buf,['  DEPEND_0          = time_tags__' DATASET_ID '\n']);
	buf = sprintf('%s%s',buf,['END_VARIABLE        = ' dsc.name{j} '__' DATASET_ID '\n!\n']);
end

buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,'!                       Data                          !\n');
buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,'DATA_UNTIL = EOF\n');

[fid,msg] = fopen([file_name ext_s],'w');
if fid < 0
	irf_log('save',['problem opening CEF file: ' msg])
	status = 1;
	return
end

sta = fprintf(fid,buf);
fclose(fid);
if sta<=0, irf_log('save','problem writing CEF header'), status = 1; return, end

if ~isempty(data)
	n_col = size(data,2) -1; % number of data columns - time
	for j=1:n_col
		ii = find(isnan(data(:,j+1)));
		if ~isempty(ii), data(ii,j+1) = FILL_VAL; end
	end
	
	s = cefprint_mx([file_name ext_s],data);
	if s~=0
		if s==1, msg = 'problem writing CEF data';
		elseif s==2, msg = 'problem compressing CEF';
		else msg = 'unknown error';
		end
		irf_log('save',msg)
		status = 1;
		return
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obuf = pmeta(buf,m_s,s,cl_id)
% Print META

obuf = sprintf('%s%s',buf,['START_META     =   ' m_s '\n']);
if iscell(s)
	for j=1:length(s)
		if isnumeric(s{j}), q = ''; ss = num2str(s{j}); else q = '"'; ss = s{j};end
		obuf = sprintf('%s%s',obuf,['   ENTRY       =   ' q ss q '\n']); 
	end
else
	if nargin==4, ss = irf_ssub(s,cl_id); 
    else ss = s;
	end
	if isnumeric(ss), q = ''; ss = num2str(ss); else q = '"'; end
	obuf = sprintf('%s%s',obuf,['   ENTRY       =   ' q ss q '\n']);
end
obuf = sprintf('%s%s',obuf,['END_META       =   ' m_s '\n']);
