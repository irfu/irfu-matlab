function status = caa_export(lev,caa_vs,cl_id,QUALITY,DATA_VERSION,sp,st,dt)
%CAA_EXPORT  export data to CAA CEF files
%
% STATUS = caa_export(data_level,caa_vs,cl_id,QUALITY,DATA_VERSION,sp,st,dt)
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
EFW_DATASET_VERSION = '3';

if nargin<8, error('time interval needed'); end    % Now REQUIRED, for caa_get ! (ML)
% The above line nullifies these:
%if nargin<6, sp='.'; end
%if nargin<5, DATA_VERSION = '02'; end
%if nargin<4, QUALITY = 3; end % Good for publication, subject to PI approval
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
			v_size = 3;    % Previously 1. (ML)
			sfit_probe = caa_sfit_probe(cl_id);
		elseif lev==3
			sfit_probe = caa_sfit_probe(cl_id);
			vs = irf_ssub('diEs?p!',cl_id,sfit_probe);
			irf_log('proc',sprintf('using p%d',sfit_probe))
			v_size = 4;    % Previously 2. (ML)
		else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'DER'
	   %keyboard
	   %vs = 'Dadc?p!';
	   vs = sprintf('Dadc%.0fp?', cl_id);
	   probe_pairs = [12, 32, 34];
%		sfit_probe = caa_sfit_probe(cl_id);
%		vs = irf_ssub('Dadc?p!',cl_id,sfit_probe);
%		irf_log('proc',sprintf('using p%d',sfit_probe))    % TODO: Change printing!
		v_size = 2;
	otherwise
		error('unknown variable')
	end
end
keyboard
% Load data
if strcmp(caa_vs, 'DER')
%   [ok, data] = c_load(vs,cl_id,'res',probe_pairs);   % Try loading data for all probe pairs.
   for k = 1:length(probe_pairs)       % Try loading data for all probe pairs.
      [data{k} ok(k)] = caa_get(st, dt, cl_id, irf_ssub(vs, probe_pairs(k)));
   end
   probe_pairs = probe_pairs(logical(ok));   % Keep list of probe pairs actually loaded.
   vs = irf_ssub(vs, probe_pairs(1));
else
%   [ok,data] = c_load(vs);
   [data, ok] = caa_get(st, dt, cl_id, vs);
%   keyboard
end
if all(~ok) || isempty(data)
	irf_log('load', ['No ' vs])
	cd(old_pwd)
	return
end
d_info = []; ok = 0;
try
%	[ok, d_info] = c_load([vs '_info'],'var');
   [d_info, ok] = caa_get(st, dt, cl_id, [vs '_info'], 'load_args', 'var');
catch
   irf_log('load', ['No ' vs '_info'])
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
			%keyboard
   if strcmp(caa_vs, 'DER')
      data1 = irf_tlim([data{1:2}], t_int);
      data2 = irf_tlim(data{3}, t_int);
      if isempty(data1) && isempty(data2)
		   irf_log('save', 'Saving empty subinterval')
	   end
   else
%      keyboard
	   data = irf_tlim(data,t_int);  % NOTE: Superfluous when using caa_get above. (ML)
	   if isempty(data)
		   irf_log('save', 'Saving empty subinterval')
	   end 
	end
else
	[iso_ts,dtint] = caa_read_interval;
	if strcmp(caa_vs, 'DER')
	   disp('Check data for DER case when no time interval given!')   % TODO: Check this!
	   keyboard
	end
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
	
	% Remove Ez, which is zero
	if lev == 2 || lev == 3
	   data = data(:, [1:3 5:end]);     % Remove column 4 (Ez data)
   else
      irf_log('warn', 'Ez not removed (data level not 2 or 3)!')
	end
	
	% Extend data array to accept bitmask and quality flag (2 columns at the end)
	data = [data zeros(size(data, 1), 2)];
   data(:, end) = QUALITY;    % Default quality column to best quality, i.e. good data/no problems.
   quality_column = size(data, 2);
   bitmask_column = quality_column - 1;
%	keyboard
	% Identify and flag problem areas in data with bitmask and quality factor:
	data = caa_identify_problems(data, lev, dsc.sen, cl_id, bitmask_column, quality_column);
	
	% Extend variable description to include the new columns bitmask and quality:
	dsc.cs = {dsc.cs{:}, 'na', 'na'};
	dsc.rep = {dsc.rep{:}, '', ''};
	dsc.units =  {dsc.units{:}, 'unitless', 'unitless'};
	dsc.si_conv = {dsc.si_conv{:}, '', ''};
	dsc.size = [dsc.size, 1, 1];
	dsc.tensor_order = [dsc.tensor_order, 0, 0];
	dsc.name = {dsc.name{:}, 'E_bitmask', 'E_quality'};
	dsc.labels = {dsc.labels{:}, 'Bitmask', 'Quality'};
	dsc.label_1 = {dsc.label_1{:}, '', ''};
	dsc.col_labels = {dsc.col_labels{:}, '', ''};
	dsc.rep_1 = {dsc.rep_1{:}, '', ''};
	dsc.field_name = {dsc.field_name{:}, ...
		'Electric field measurement quality bitmask',...
		'Electric field measurement quality flag (9=best)'};
	dsc.ptype = {dsc.ptype{:}, 'Support_Data', 'Support_Data'};
	dsc.valtype = {dsc.valtype{:}, 'INT', 'INT'};
	dsc.sigdig = [dsc.sigdig, 5, 1];
	dsc.ent = {dsc.ent{:}, 'Electric_Field', 'Electric_Field'};
	dsc.prop = {dsc.prop{:}, 'Status', 'Status'};
	dsc.fluc = {dsc.fluc{:}, '', ''};
	
	
	% Correct offsets
	if ~isempty(data)
		
		if lev==3
			% Delta offsets: remove automatic and apply CAA
			Del_caa = c_efw_delta_off(data(1,1),cl_id);
			if ~isempty(Del_caa)
%				[ok,Delauto] = c_load('D?p12p34',cl_id);
            [Delauto, ok] = caa_get(st, dt, cl_id, 'D?p12p34');
				if ~ok || isempty(Delauto)
					irf_log('load',irf_ssub('Cannot load/empty D?p12p34',cl_id))
				else
					data = caa_corof_delta(data,sfit_probe,Delauto,'undo');
					data = caa_corof_delta(data,sfit_probe,Del_caa,'apply');
				end
			end
		end
		
		% Dsi offsets
		dsiof = c_ctl(cl_id,'dsiof');
		if isempty(dsiof)
%			[ok,Ps,msg] = c_load('Ps?',cl_id);
         [Ps, ok] = caa_get(st, dt, cl_id, 'Ps?');
%			if ~ok, irf_log('load',msg), end
         if ~ok, irf_log('load',irf_ssub('Cannot load/empty Ps?',cl_id)), end
			if caa_is_sh_interval
				[dsiof_def, dam_def] = c_efw_dsi_off(t_int(1),cl_id,[]);
			else
				[dsiof_def, dam_def] = c_efw_dsi_off(t_int(1),cl_id,Ps);
			end

%			[ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end
         [Ddsi, ok1] = caa_get(st, dt, cl_id, 'Ddsi?'); if ~ok1, Ddsi = dsiof_def; end
%			[ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end
			[Damp, ok2] = caa_get(st, dt, cl_id, 'Damp?'); if ~ok2, Damp = dam_def; end
			
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
		
		%if length(Ddsi) == 1
		%	dsc.com = sprintf('ISR2 offsets: dEx=%1.2f dEy=%1.2f dAmp=%1.2f',...
		%		real(Ddsi(1)),imag(Ddsi(1)),Damp);
		%else
		%	dsc.com = 'ISR2 offsets';
		%	for in = 1:size(Ddsi,1)
		%		dsc.com = [dsc.com sprintf(' %s: dEx=%1.2f dEy=%1.2f,',...
		%			epoch2iso(Ddsi(in,1),1),real(Ddsi(in,2)),imag(Ddsi(in,2)))];
		%	end
		%	dsc.com = [dsc.com sprintf(' dAmp=%1.2f',Damp)];
		%end
		clear Ddsi Damp
		
		dsc.com = 'For offsets see DER dataset for this interval';
		dsc.com = [dsc.com '. Probes: ' dsc.sen];
		irf_log('calb',dsc.com)
		
	end
	
%	[ok,Del] = c_load('D?p12p34',cl_id); if ~ok, Del = [0 0]; end
   [Del, ok] = caa_get(st, dt, cl_id, 'D?p12p34'); if ~ok, Del = [0 0]; end
%	if ~isreal(Del)
%		Del = imag(Del);
%		% offset is applied to p34
%		if strcmp(dsc.sen,'12') || strcmp(dsc.sen,'32'), Del = [0 0]; end
%		dsc.com = sprintf('%s. p%s offset (ISR2): dEx=%1.2f dEy=%1.2f',...
%			dsc.com, '34', Del(1), Del(2));
%	else
%		% offset is applied to p12/32
%		if strcmp(dsc.sen,'34'), Del = [0 0]; end
%		dsc.com = sprintf('%s. p%s offset (ISR2): dEx=%1.2f dEy=%1.2f',...
%			dsc.com, dsc.sen(1:2), Del(1), Del(2));
%	end
	
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
	
% Combine ADC offsets from two probe pairs into one dataset:
elseif strcmp(caa_vs, 'DER')

   if length(probe_pairs) == 1   % Most likely p1 broken, so no 'Dadc?p12' product!
      single_pair_data = [data1 data2];
      start_time = min( single_pair_data(:,1) );
   else
%      start_time = min( min(data1(:,1)), min(data2(:,1)) );
      start_time = min( min( [data1(:,1) data2(:,1)] ) );
   end
   timestamp = start_time:4:t_int(2);
   data_out = zeros(length(timestamp), 3) * NaN;
   data_out(:, 1) = timestamp;

   if find(probe_pairs == 12 | probe_pairs == 32)
      [ind1, ind2] = irf_find_comm_idx(data_out, data1);
      data_out(ind1, 2) = data1(ind2, 2);
   end
   
   if find(probe_pairs == 34)
      [ind1, ind2] = irf_find_comm_idx(data_out, data2);
      data_out(ind1, 3) = data2(ind2, 2);
   end
   
   data_orig = data;
   clear data;
   data = data_out;
   
   % Extend description to cover data record for two probe pairs:
   if length(probe_pairs) > 1
      dsc.com = sprintf('Probe pairs used are p%i and p%i', probe_pairs);
   else
      dsc.com = sprintf('Probe pair used is p%i', probe_pairs);
   end
      dsc.valtype = {dsc.valtype{:}, 'FLOAT'};
      dsc.sigdig = [dsc.sigdig 6];
      dsc.size = [dsc.size 1];
   
   clear start_time timestamp ind1 ind2 data_out
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
buf = sprintf('%s%s',buf,irf_ssub('include = "C?_CH_EFW_INST.ceh"\n',cl_id));
buf = sprintf('%s%s',buf,irf_ssub('include = "C?_CH_EFW_L!_$.ceh"\n', ...
               cl_id, lev, caa_vs)); % Change to 'E', 'P', etc.!
buf = pmeta(buf,'FILE_TYPE','cef');
buf = pmeta(buf,'DATASET_VERSION',EFW_DATASET_VERSION);
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
buf = pmeta(buf, 'FILE_CAVEATS', dsc.com);


buf = sprintf('%s%s',buf,'!\n');
buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,'!                       Data                          !\n');
buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,'DATA_UNTIL = "END_OF_DATA"\n');

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
	
	% build up formatting string
    format=repmat('%8.3f',n_col,1);
    i=1;
    for j=1:v_size
        if strcmp(dsc.valtype{j},'FLOAT')==0
            fcode=['%' num2str(dsc.sigdig(j)+1,'%1.0f')  '.0f'];
            format(i:i+dsc.size(j)-1,:) = repmat(fcode,dsc.size(j),1);
        end
        i=i+dsc.size(j);
    end
    format=ctranspose(format);
	
   s = cefprint_mx([file_name ext_s],data, format);
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
