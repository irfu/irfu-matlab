function caa_export(lev,caa_vs,cl_id,QUALITY,DATA_VERSION,sp,st,dt)
%CAA_EXPORT  export data to CAA CEF files
%
% caa_export(data_level,caa_vs,cl_id,[QUALITY,DATA_VERSION,sp,st,dt])
%
% CEF file is written to a current directory. Data is supposed to be also 
% there, if SP is not specified.
%
% See also C_EXPORT_ASCII
%
% $Id$

% Copyright 2004-2006 Yuri Khotyaintsev

% This must be changed when we do any major changes to our processing software
EFW_DATASET_VERSION = '1';

if nargin<8, st = []; dt=[]; end
if nargin<6, sp='.'; end
if nargin<5, DATA_VERSION = '01'; end
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
		elseif lev==3
			sfit_probe = caa_sfit_probe(cl_id);
			vs = irf_ssub('diEs?p!',cl_id,sfit_probe);
			irf_log('proc',sprintf('using p%d',sfit_probe))
			v_size = 2;
		else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'EF'
		vs = irf_ssub('diEF?p1234',cl_id);
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
end

if ~ok || isempty(d_info), dsc = c_desc(vs);
else dsc = c_desc(vs,d_info);
end

if lev==3
	TIME_RESOLUTION = 4;
elseif (lev==1 && ~isempty(regexp(caa_vs,'^P(1|2|3|4)?$'))) || ...
		(lev==2 && strcmp(caa_vs,'P'))
	TIME_RESOLUTION = 1/5;
elseif (lev==1 && ~isempty(regexp(caa_vs,'^P(12|32|34)?$'))) || ...
		(lev==2 && (strcmp(caa_vs,'E') || strcmp(caa_vs,'EF')))
	fs = c_efw_fsample(data,'hx');
	if ~fs, error('cannot determine time resolution'), end
	TIME_RESOLUTION = 1/fs;
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
	% We check if this full res E is from coming from two probe pairs
	if lev==2 && ~(strcmp(dsc.sen,'1234') || strcmp(dsc.sen,'3234')) && QUALITY>1
		irf_log('save','This is not a full E, setting QUALITY=1!')
		QUALITY = 1;
	end
	
	if ~isempty(data)
		dsiof = c_ctl(cl_id,'dsiof');
		if isempty(dsiof), dsiof = [1+0i 1]; end
		[ok,Ddsi] = c_load('Ddsi?',cl_id); if ~ok, Ddsi = dsiof(1); end
		[ok,Damp] = c_load('Damp?',cl_id); if ~ok, Damp = dsiof(2); end
		clear dsiof
		
		data = caa_corof_dsi(data,Ddsi,Damp);
		dsc.com = sprintf('ISR2 offsets: dEx=%1.2f dEy=%1.2f dAmp=%1.2f',...
			real(Ddsi(1)),imag(Ddsi(1)),Damp);
		irf_log('calb',dsc.com)
		dsc.com = [dsc.com '. Probes: ' dsc.sen];
		
		% Remove Ez, which is zero
		if lev==3, data = data(:,[1:3 5]);
        else data = data(:,1:3);
		end
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
	
elseif lev==1 && ~isempty(regexp(caa_vs,'^P(12|32|34)?$'))
	if ~isempty(data)
		% convert mV/m back to V
		if id==32, data(:,2) = data(:,2)*.0622;
        else data(:,2) = data(:,2)*.088;
		end
		id = str2double(caa_vs(2:end));
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
fid = fopen([file_name ext_s],'w');

fprintf(fid,'!-------------------- CEF ASCII FILE --------------------|\n');
nnow = now;
fprintf(fid,['! created on ' datestr(nnow) '\n']);
fprintf(fid,'!--------------------------------------------------------|\n');
fprintf(fid,['FILE_NAME = "' file_name ext_s '"\n']);
fprintf(fid,'FILE_FORMAT_VERSION = "CEF-2.0"\n');
fprintf(fid,['END_OF_RECORD_MARKER = "' EOR_MARKER '"\n']);
fprintf(fid,'include = "CL_CH_MISSION.ceh"\n');
fprintf(fid,irf_ssub('include = "C?_CH_OBS.ceh"\n',cl_id));
fprintf(fid,'include = "CL_CH_EFW_EXP.ceh"\n');
pmeta(fid,'FILE_TYPE','cef')
pmeta(fid,'DATA_TYPE','CP')
pmeta(fid,'INSTRUMENT_NAME','EFW?',cl_id)
pmeta(fid,'INSTRUMENT_DESCRIPTION','EFW Experiment on Cluster C?',cl_id)
pmeta(fid,'INSTRUMENT_CAVEATS','*C?_CQ_EFW_CAVEATS',cl_id)
pmeta(fid,'DATASET_ID',DATASET_ID,cl_id)
pmeta(fid,'DATASET_TITLE',dsc.field_name{1})
pmeta(fid,'DATASET_DESCRIPTION',...
	{'This dataset contains measurements of the', ...
	[DATASET_DESCRIPTION_PREFIX dsc.field_name{1}],... 
	irf_ssub('from the EFW experiment on the Cluster C? spacecraft',cl_id)})
pmeta(fid,'DATASET_VERSION',EFW_DATASET_VERSION)
pmeta(fid,'TIME_RESOLUTION',TIME_RESOLUTION)
pmeta(fid,'MIN_TIME_RESOLUTION',TIME_RESOLUTION)
pmeta(fid,'MAX_TIME_RESOLUTION',TIME_RESOLUTION)
pmeta(fid,'PROCESSING_LEVEL',PROCESSING_LEVEL)
pmeta(fid,'DATASET_CAVEATS',['*C?_CQ_EFW_' caa_vs],cl_id)
pmeta(fid,'LOGICAL_FILE_ID',file_name)
pmeta(fid,'VERSION_NUMBER',DATA_VERSION)
fprintf(fid,'START_META     =   FILE_TIME_SPAN\n');
fprintf(fid,'   VALUE_TYPE  =   ISO_TIME_RANGE\n');
fprintf(fid,...
['   ENTRY       =   ' epoch2iso(t_int(1)) '/' epoch2iso(t_int(2)) '\n']);
fprintf(fid,'END_META       =   FILE_TIME_SPAN\n');
fprintf(fid,'START_META     =   GENERATION_DATE\n');
fprintf(fid,'   VALUE_TYPE  =   ISO_TIME\n');
fprintf(fid,['   ENTRY       =   ' epoch2iso(date2epoch(nnow)) '\n']);
fprintf(fid,'END_META       =   GENERATION_DATE\n');

fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fid,'!                   Variables                         !\n');
fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fid,['START_VARIABLE    = time_tags__' DATASET_ID '\n']);
fprintf(fid,'  VALUE_TYPE      = ISO_TIME\n');
fprintf(fid,['  DELTA_PLUS      = ' num2str(TIME_RESOLUTION/2) '\n']);
fprintf(fid,['  DELTA_MINUS     = ' num2str(TIME_RESOLUTION/2) '\n']);
fprintf(fid, '  FILLVAL         = 9999-12-31T23:59:59Z\n');
fprintf(fid,'  LABLAXIS        = "UT"\n');
fprintf(fid,'  FIELDNAM        = "Universal Time"\n');
fprintf(fid,['END_VARIABLE      = time_tags__' DATASET_ID '\n!\n']);

for j=1:v_size
	fprintf(fid,['START_VARIABLE      = ' dsc.name{j} '__' DATASET_ID '\n']);
	fprintf(fid,'  PARAMETER_TYPE    = "Data"\n');
	fprintf(fid,'  SIZES             = %d\n',dsc.size(j));
	fprintf(fid,'  VALUE_TYPE        = FLOAT\n');
	fprintf(fid,['  ENTITY            = "' dsc.ent{j} '"\n']);
	fprintf(fid,['  PROPERTY          = "' dsc.prop{j} '"\n']);
	if ~isempty(dsc.fluc{j})
		fprintf(fid,['  FLUCTUATIONS      = "' dsc.fluc{j} '"\n']);
	end
	fprintf(fid,['  CATDESC           = "' dsc.field_name{j} '"\n']);
	fprintf(fid,['  FIELDNAM          = "' dsc.field_name{j} '"\n']);
	if ~isempty(dsc.si_conv{j})
		fprintf(fid,['  SI_CONVERSION     = "' dsc.si_conv{j} '"\n']);
	else
		fprintf(fid,['  SI_CONVERSION     = "1>' dsc.units{j} '"\n']);
	end
	fprintf(fid,['  UNITS             = "' dsc.units{j} '"\n']);
	fprintf(fid,['  FILLVAL           = ' num2str(FILL_VAL,'%8.3f') '\n']);
	fprintf(fid,['  QUALITY           = ' num2str(QUALITY) '\n']);
	fprintf(fid,'  SIGNIFICANT_DIGITS= 6 \n');
	if ~isempty(dsc.com) && j==1
		fprintf(fid,['  PARAMETER_CAVEATS = "' dsc.com '"\n']);
	end
	if ~strcmp(dsc.cs{j},'na')
		fprintf(fid,['  COORDINATE_SYSTEM = "' dsc.cs{j} '"\n']); 
	end
	if isfield(dsc,'frv')
		if ~isempty(dsc.frv{j})
			fprintf(fid,['  FRAME_VELOCITY    = "' dsc.frv{j} '"\n']);
		end
	end
	if dsc.size(j) > 1
		fprintf(fid,['  REPRESENTATION_1  = ' dsc.rep_1{j} '\n']);
		fprintf(fid,['  LABEL_1           = ' dsc.label_1{j} '\n']);
	end
	fprintf(fid,['  LABLAXIS          = "' dsc.labels{j} '"\n']);
	fprintf(fid,['  DEPEND_0          = time_tags__' DATASET_ID '\n']);
	fprintf(fid,['END_VARIABLE        = ' dsc.name{j} '__' DATASET_ID '\n!\n']);
end

fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fid,'!                       Data                          !\n');
fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fid,'DATA_UNTIL = EOF\n');

fclose(fid);

if ~isempty(data)
	n_col = size(data,2) -1; % number of data columns - time
	for j=1:n_col
		ii = find(isnan(data(:,j+1)));
		if ~isempty(ii), data(ii,j+1) = FILL_VAL; end
	end
	
	s = cefprint_mx([file_name ext_s],data);
	if s~=0, irf_log('save','problem writing to CEF file'), return, end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pmeta(fid,m_s,s,cl_id)
% Print META

fprintf(fid,['START_META     =   ' m_s '\n']);
if iscell(s)
	for j=1:length(s)
		if isnumeric(s{j}), q = ''; ss = num2str(s{j}); else q = '"'; ss = s{j};end
		fprintf(fid,['   ENTRY       =   ' q ss q '\n']); 
	end
else
	if nargin==4, ss = irf_ssub(s,cl_id); 
    else ss = s;
	end
	if isnumeric(ss), q = ''; ss = num2str(ss); else q = '"'; end
	fprintf(fid,['   ENTRY       =   ' q ss q '\n']);
end
fprintf(fid,['END_META       =   ' m_s '\n']);
