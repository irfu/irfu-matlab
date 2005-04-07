function caa_export(lev,caa_vs,cl_id,DATA_VERSION,sp)
%CAA_EXPORT export data to CAA CEF files
% caa_export(lev,caa_vs,cl_id,DATA_VERSION,sp)
%
% See also c_export_ascii
%
% $Id$

% Copyright 2004,2005 Yuri Khotyaintsev

% This must be changed when we do any major changes to our processing software
EFW_DATASET_VERSION = '1';
QUALITY = '3'; % Good for publication, subject to PI approval

if nargin<5, sp='.'; end
if nargin<4, DATA_VERSION = '01'; end
if cl_id<=0 | cl_id>4, error('CL_ID must be 1..4'), end
if lev<1 | lev>3, error('LEV must be 1,2 or 3'), end

DATASET_DESCRIPTION_PREFIX = '';

if lev==1
	if regexp(caa_vs,'^P(1|2|3|4|12|32|34)?$')
		id = str2num(caa_vs(2:end));
		if id <=4, vs = irf_ssub('P10Hz?p!',cl_id,id);
		else, vs = irf_ssub('wE?p!',cl_id,id);
		end
		v_size = 1;
	else, error('Must be P(1|2|3|4|12|32|34)')
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
			vs = irf_ssub('diEs?p34',cl_id);
			v_size = 2;
		else
			disp('not implemented'), return
		end
	case 'EF'
		disp('not implemented'), return
	otherwise
		error('unknown variable')
	end
end

EOR_MARKER = '$';
FILL_VAL = -1.0E9;
PROCESSING_LEVEL='Calibrated';

old_pwd = pwd;
cd(sp)

% Load data
data = c_load(vs,'var');
if isempty(data)
	irf_log('load', ['No ' vs])
	cd(old_pwd)
	return
end
d_info = [];
try
	d_info = c_load([vs '_info'],'var');
end

if isempty(d_info), dsc = c_desc(vs);
else, dsc = c_desc(vs,d_info);
end

dt = data(2,1) - data(1,1);
ddt = .05;
dt_spin = 4;
dt_lx = 1/5;
dt_nm = 1/25;
dt_bm1 = 1/450;

if dt>(1-ddt)*dt_spin & dt<(1+ddt)*dt_spin
	TIME_RESOLUTION = dt_spin;
elseif dt>(1-ddt)*dt_lx & dt<(1+ddt)*dt_lx
	TIME_RESOLUTION = dt_lx;
elseif dt>(1-ddt)*dt_nm & dt<(1+ddt)*dt_nm
	TIME_RESOLUTION = dt_nm;
elseif dt>(1-ddt)*dt_bm1 & dt<(1+ddt)*dt_bm1
	TIME_RESOLUTION = dt_bm1;
else,
	error('cannot determine time resolution')
end

% Do magic on E-field
if strcmp(caa_vs,'E')
	% We check if this full res E is from coming from two probe pairs
	if lev==2 & ~(strcmp(dsc.sen,'1234') | strcmp(dsc.sen,'3234'))
		c_log('save','This is not a full E, we skip it.')
		return
	end
	
	Dxy_s =  irf_ssub('Ddsi?',cl_id);
	Dx_s =  irf_ssub('real(Ddsi?)',cl_id);
	Dy_s =  irf_ssub('imag(Ddsi?)',cl_id);
	Da_s =  irf_ssub('Damp?',cl_id);

	eval(['load mEDSI ' Dxy_s ' ' Da_s])
	if exist(Dxy_s,'var'), eval(['Dx=real(' Dxy_s ');Dy=imag(' Dxy_s ');'])
	else, irf_log('calb','using Dx=1,Dy=0'), Dx = 1; Dy=0;
	end
	if exist(Da_s,'var'), eval(['Da=' Da_s ';'])
	else, disp('using Da=1'), Da = 1;
	end

	data = caa_corof_dsi(data,Dx,Dy,Da);
	dsc.com = sprintf('ISR2 offsets: dEx=%1.2f dEy=%1.2f dAmp=%1.2f',Dx,Dy,Da);
	dsc.com = [dsc.com '. Probes: ' dsc.sen];
	
	%dsc.fro = {'Coordinate_System'};
	dsc.frv = {'Observatory'};
	if v_size>1
		for j=2:v_size
			%dsc.fro = [dsc.fro {''}]; 
			dsc.frv = [dsc.frv {''}];
		end
	end
	
	% Remove Ez, which is zero
	if lev==3, data = data(:,[1:3 5]);
	else, data = data(:,1:3);
	end
elseif lev==1 & regexp(caa_vs,'^P(12|32|34)?$')
	% convert mV/m back to V
	id = str2num(caa_vs(2:end));
	
	dsc.units = {'V'};
	
	if id==32, data(:,2) = data(:,2)*.0622;
	else, data(:,2) = data(:,2)*.088;
	end
	
	dsc.cs = {'na'};
 	dsc.units =  {'V'};
	dsc.si_conv = {''};
	dsc.field_name = {['Potential difference measured between probes '...
		dsc.sen(1) ' and ' dsc.sen(2)]};
	dsc.ent = {'Instrument'};
	dsc.prop = {'Probe_Potential'};
	dsc.fluc = {'Waveform'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write to file
ext_s = '.cef';
% We have special names for CAA
DATASET_ID = irf_ssub(['C?_CP_EFW_L' num2str(lev) '_' caa_vs],cl_id);
file_name = ...
	[DATASET_ID '__' irf_fname(data([1 end],1),2) '_V' DATA_VERSION];
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
%pmeta(fid,'OBSERVATORY','Cluster-?',cl_id)
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
pmeta(fid,'TIME_RESOLUTION',num2str(TIME_RESOLUTION))
pmeta(fid,'MIN_TIME_RESOLUTION',num2str(TIME_RESOLUTION))
pmeta(fid,'MAX_TIME_RESOLUTION',num2str(TIME_RESOLUTION))
pmeta(fid,'PROCESSING_LEVEL',PROCESSING_LEVEL)
pmeta(fid,'DATASET_CAVEATS',['*C?_CQ_EFW_' caa_vs],cl_id)
pmeta(fid,'LOGICAL_FILE_ID',file_name)
pmeta(fid,'VERSION_NUMBER',DATA_VERSION)
fprintf(fid,'START_META     =   FILE_TIME_SPAN\n');
fprintf(fid,'   VALUE_TYPE  =   ISO_TIME_RANGE\n');
fprintf(fid,...
['   ENTRY       =   "' epoch2iso(data(1,1)) '/' epoch2iso(data(end,1)) '"\n']);
fprintf(fid,'END_META       =   FILE_TIME_SPAN\n');
fprintf(fid,'START_META     =   GENERATION_DATE\n');
fprintf(fid,'   VALUE_TYPE  =   ISO_TIME\n');
fprintf(fid,['   ENTRY       =   "' epoch2iso(date2epoch(nnow)) '"\n']);
fprintf(fid,'END_META       =   GENERATION_DATE\n');

fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fid,'!                   Variables                         !\n');
fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fid,['START_VARIABLE    = time_tags__' DATASET_ID '\n']);
fprintf(fid,'  VALUE_TYPE      = ISO_TIME\n');
fprintf(fid,['  DELTA_PLUS      = ' num2str(TIME_RESOLUTION/2) '\n']);
fprintf(fid,['  DELTA_MINUS     = ' num2str(TIME_RESOLUTION/2) '\n']);
fprintf(fid,['  FILLVAL         = "9999-12-31T23:59:59Z"\n']);
fprintf(fid,'  LABLAXIS        = "UT"\n');
fprintf(fid,'  FIELDNAM        = "Universal Time"\n');
fprintf(fid,['END_VARIABLE      = time_tags__' DATASET_ID '\n!\n']);

for j=1:v_size
	fprintf(fid,['START_VARIABLE      = ' dsc.name{j} '__' DATASET_ID '\n']);
	fprintf(fid,'  PARAMETER_TYPE    = "Data"\n');
	fprintf(fid,'  SIZES             = %d\n',dsc.size(j));
	fprintf(fid,'  VALUE_TYPE        = FLOAT\n');
	fprintf(fid,'  MEASUREMENT_TYPE  = "Electric_Field"\n');
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
	fprintf(fid,['  FILLVAL           = "' num2str(FILL_VAL,'%8.3f') '"\n']);
	fprintf(fid,['  QUALITY           = "' QUALITY '"\n']);
	fprintf(fid,'  SIGNIFICANT_DIGITS= 6 \n');
	if ~isempty(dsc.com) & j==1
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

time_s = epoch2iso(data(:,1));
sz = size(data);
n_data = sz(2) - 1; % number of data columns - time
mask = '';
for j=1:n_data
	mask = [mask ', %8.3f'];
	ii = find(isnan(data(:,j+1)));
	if ~isempty(ii), data(ii,j+1) = FILL_VAL; end
end

for j=1:length(data(:,1))
	fprintf(fid,time_s(j,:));
	fprintf(fid,[mask ' ' EOR_MARKER '\n'],data(j,2:end)');
end 


fclose(fid);
cd(old_pwd)

function pmeta(fid,m_s,s,cl_id)
% Print META

fprintf(fid,['START_META     =   ' m_s '\n']);
if iscell(s)
	for j=1:length(s), fprintf(fid,['   ENTRY       =   "' s{j} '"\n']); end
else
	if nargin==4, ss = irf_ssub(s,cl_id); 
	else, ss = s;
	end
	fprintf(fid,['   ENTRY       =   "' ss '"\n']);
end
fprintf(fid,['END_META       =   ' m_s '\n']);
