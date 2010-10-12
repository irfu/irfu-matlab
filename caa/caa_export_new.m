function status = caa_export_new(lev,caa_vs,cl_id,QUALITY,DATA_VERSION,sp,st,dt)
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

%if nargin<8, st = []; dt=[]; end
if nargin<8, error('time interval needed'); end    % Now REQUIRED, for caa_get ! (ML)
% The above line nullifies these:
%if nargin<6, sp='.'; end
%if nargin<5, DATA_VERSION = '00'; end
%if nargin<4, QUALITY = 3; end % Good for publication, subject to PI approval
if cl_id<=0 || cl_id>4, error('CL_ID must be 1..4'), end
if lev<1 || lev>3, error('LEV must be 1,2 or 3'), end

DATASET_DESCRIPTION_PREFIX = '';
EOR_MARKER = '$';
FILL_VAL = -1.0E9;
PROCESSING_LEVEL='Calibrated';
DELIVERY_TO_CAA = 1;    % Changes file name format to new CAA daily-file-format

old_pwd = pwd;
%cd(sp)
dirs = caa_get_subdirs(st, dt, cl_id);
%if isempty(dirs), disp(['Invalid (empty?) dir: ' sp]), cd(old_pwd); status = 1; return, end

if DELIVERY_TO_CAA   % Check that files are midnight-to-midnight
   st_temp = fromepoch(st);
   if dt ~= 24*3600 || st_temp(4) ~= 00
      error('Incorrect format for daily files. Check parameters!')
   end
end

result = [];
result_com = {};


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
%			sfit_probe = caa_sfit_probe(cl_id);
		elseif lev==3
%			sfit_probe = caa_sfit_probe(cl_id);
%			vs = irf_ssub('diEs?p!',cl_id,sfit_probe);
%			irf_log('proc',sprintf('using p%d',sfit_probe))
			v_size = 4;    % Previously 2. (ML)
		else
			disp('not implemented'), cd(old_pwd), return
        end
	case 'HK'
		if lev==2
			vs = irf_ssub('HK?',cl_id);
			v_size = 12;
		else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'DER'
	   %keyboard
	   %vs = 'Dadc?p!';
	   vs_DER = sprintf('Dadc%.0fp?', cl_id);
	   probe_pair_list = [12, 32, 34];
%		sfit_probe = caa_sfit_probe(cl_id);
%		vs = irf_ssub('Dadc?p!',cl_id,sfit_probe);
%		irf_log('proc',sprintf('using p%d',sfit_probe))    % TODO: Change printing!
		v_size = 2;
	case 'SFIT'
		if lev==3
            % Fake for c_desc only. No data variable in .mat files
			vs = irf_ssub('SFIT?',cl_id);
			v_size = 4;
		else
			disp('not implemented'), cd(old_pwd), return
		end
	otherwise
		error('unknown variable')
	end
end


for dd = 1:length(dirs)
   d = dirs{dd};
%   cd([sp '/' d]);
   cd(d);
   
   % Set up spin fit related information.
   % Probe pair used can vary for each subinterval!
   if strcmp(caa_vs, 'E')
      sfit_probe = caa_sfit_probe(cl_id);
      irf_log('proc',sprintf('using p%d',sfit_probe))
      if lev == 3
         vs = irf_ssub('diEs?p!',cl_id,sfit_probe);
      end
   end

   %   keyboard
   % Load data
   if strcmp(caa_vs, 'DER')
      vs = vs_DER;
      [ok, data] = c_load(vs,cl_id,'res',probe_pair_list);   % Try loading data for all probe pairs.
   %   for k = 1:length(probe_pair_list)       % Try loading data for all probe pairs.
   %      [data{k} ok(k)] = caa_get(st, dt, cl_id, irf_ssub(vs, probe_pair_list(k)));
   %   end
      probe_pairs = probe_pair_list(logical(ok));   % Keep list of probe pairs actually loaded.
      if numel(probe_pairs) > 0
         vs = irf_ssub(vs, probe_pairs(1));
      end
   elseif strcmp(caa_vs, 'SFIT')
      [ok,spf34,msg] = c_load('diEs?p34',cl_id);
      if ~ok || isempty(spf34)
		 irf_log('load',msg)
		 data = []; cd(old_pwd); return
      end;
      nanfill = 0;
      pnosfit = 12;
      ret=whos('-file','./mEDSI.mat',irf_ssub('diEs?p!',cl_id,pnosfit));
      if isempty(ret)
         pnosfit = 32;
      end
      [ok,spfD,msg] = c_load(irf_ssub('diEs?p!',cl_id,pnosfit));
      if ~ok || isempty(spfD)
         nanfill = 1; % No P12/32 data
         irf_log('load',irf_ssub('Fillvalue used. No diEs?p12 or diEs?p32 data',cl_id) )
         ok = 1; % yes it's all good
      else
         [ok,Del,msg] = c_load('D?p12p34',cl_id);
         if ~ok || isempty(Del)
            irf_log('load',msg)
            data = []; cd(old_pwd); return
         end;
         if imag(Del(1)) ~= 0 || imag(Del(2)) ~= 0
            irf_log('load','Info: Imaginary delta offset.');
         end
      end;
      if nanfill
         % save time, NaN(fillval) and p34 spin-fit: B C sdev 
         data=[spf34(:,1) NaN(size(spf34,1),3) spf34(:,2:3) spf34(:,5)];
      else
         % save time, p12/32 and p34 spin-fit: B C sdev 
         data=[spf34(:,1) spfD(:,2:3)+ones(size(spfD,1),1)*Del spfD(:,5) spf34(:,2:3) spf34(:,5)];
      end
   else
      [ok,data] = c_load(vs);
   %   [data, ok] = caa_get(st, dt, cl_id, vs);
   %   keyboard
   end
   if all(~ok) || isempty(data)
   	irf_log('load', ['No ' vs]);
   	cd(old_pwd)
%   	return
      continue
   end
   d_info = []; ok = 0;
   try
   	[ok, d_info] = c_load([vs '_info'],'var');
   %   [d_info, ok] = caa_get(st, dt, cl_id, [vs '_info'], 'load_args', 'var');
   catch
      irf_log('load', ['No ' vs '_info']);
   	d_info = []; ok = 0;
   end

   if ~ok || isempty(d_info), dsc = c_desc(vs);
   else dsc = c_desc(vs,d_info);
   end
   
   if lev==3
   	TIME_RESOLUTION = 4;
   elseif (lev==1 && ~isempty(regexp(caa_vs,'^P(1|2|3|4)?$','once'))) || ...
   		(lev==2 && strcmp(caa_vs,'P'))
   	TIME_RESOLUTION = 1/5;
   elseif (lev==1 && ~isempty(regexp(caa_vs,'^P(12|32|34)?$','once'))) || ...
   		(lev==2 && strcmp(caa_vs,'E'))
   	fs = c_efw_fsample(data, 'hx', cl_id); % Use SC number for better accuracy.
   	if ~fs, error('cannot determine time resolution'), end
   	TIME_RESOLUTION = 1/fs;
   end
   
   % Make subinterval
   if ~isempty(st) && ~isempty(dt)  %  == ~(isempty(st) || isempty(dt)) == NOR
   	t_int_full = st + [0 dt];
   	
   	[iso_ts,dtint] = caa_read_interval;
   	if isempty(iso_ts)
   		t_int = data([1 end],1);
   	else
   		t_int(1) = iso2epoch(iso_ts);
   		t_int(2) = t_int(1) + dtint;
   	end
   	
   	
   	irf_log('save', sprintf('%s : %s -- %s',...
   			vs, epoch2iso(t_int(1),1), epoch2iso(t_int(2),1)))

      if strcmp(caa_vs, 'DER')
         % Limit data to both time interval given as input, and time interval read from file.
         data1 = irf_tlim([data{1:2}], t_int_full);
         data1 = irf_tlim(data1, t_int);
         data2 = irf_tlim(data{3}, t_int_full);
         data2 = irf_tlim(data2, t_int);
         if isempty(data1) && isempty(data2)
   		   irf_log('save', 'Saving empty subinterval')
   	   end
      else
         % Limit data to both time interval given as input, and time interval read from file.
   	   data = irf_tlim(data,t_int_full);  % NOTE: Superfluous when using caa_get above. (ML)
   	   data = irf_tlim(data, t_int);
   	   if isempty(data)
   		   irf_log('save', 'Saving empty subinterval')
   	   end 
   	end
   else  % isempty(st) || isempty(dt)   == ~NOR == OR
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
   if strcmp(caa_vs,'E') && ~isempty(data)
   	
   	% We export only X and Y, no need to export zeroes in Ez.
   	dsc.size(1) = 2;
   	
   	% Remove Ez, which is zero
   	data = data(:, [1:3 5:end]);     % Remove column 4 (Ez data)
   	
   	% Get info on probe pair(s) in (sub)interval.
   	if lev == 2
   	   E_info = c_load('diE?p1234_info', cl_id, 'var');  % Load info; need list of probe pairs!
   	   if isempty(E_info) || ~isfield(E_info, 'probe')
   	      error('Could not load probe pair info!')
   	   end
   	   probe_info = E_info.probe;
   	elseif lev == 3
   	   probe_info = num2str(sfit_probe);
   	end
   	
   	% Fill gap in data at start of subinterval
    %if ~isempty(data) && ~isempty(result) && ((data(1,1) - result(end,1)) > TIME_RESOLUTION)
   	%   [tmp, filldata] = caa_fill_gaps(result, data(1,1));
   	%   if ~isempty(filldata)
   	%      data = [filldata(:,1:size(data,2)); data];
   	%   end
   	%   clear tmp
   	%end
   	
   	% Fill gap in data at end of subinterval
   	%if ~isempty(data) && ((t_int(2) - data(end,1)) > TIME_RESOLUTION)
   	%   data = caa_fill_gaps(data, t_int(2));
   	%end
   	
   	% Extend data array to accept bitmask and quality flag (2 columns at the end)
   	data = [data zeros(size(data, 1), 2)];
      data(:, end) = QUALITY;    % Default quality column to best quality, i.e. good data/no problems.
      quality_column = size(data, 2);
      bitmask_column = quality_column - 1;

   	% Identify and flag problem areas in data with bitmask and quality factor:
   	data = caa_identify_problems(data, lev, probe_info, cl_id, bitmask_column, quality_column);
   	
   	% Extend variable description to include the new columns bitmask and quality:
   	dsc.size = [dsc.size, 1, 1];
   	dsc.valtype = {dsc.valtype{:}, 'INT', 'INT'};
   	dsc.sigdig = [dsc.sigdig, 5, 1];
   	
   	if length(probe_info) > 2
   	   probe_str = sprintf('Probe pairs p%s,p%s', probe_info(1:2), probe_info(3:4));
   	   if length(probe_info) > 4, probe_str = [probe_str ',p' probe_info(5:6)]; end
   	else probe_str = sprintf('Probe pair p%s', probe_info);
   	end
   	com_str = sprintf('%s/%s %s', ...
   	   epoch2iso(t_int(1),1), epoch2iso(t_int(2),1), probe_str);
   	
   	result_com{end+1} = com_str;
      irf_log('calb',com_str)
   	   
   	
   	
   	% Correct offsets
   	if ~isempty(data)
   		
   		if lev==3
   			% Delta offsets: remove automatic and apply CAA
   			Del_caa = c_efw_delta_off(data(1,1),cl_id);
   			if ~isempty(Del_caa)
   				[ok,Delauto] = c_load('D?p12p34',cl_id);
   %            [Delauto, ok] = caa_get(st, dt, cl_id, 'D?p12p34');
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
   			[ok,Ps,msg] = c_load('Ps?',cl_id);
   %         [Ps, ok] = caa_get(st, dt, cl_id, 'Ps?');
   			if ~ok, irf_log('load',msg), end
%            if ~ok, irf_log('load',irf_ssub('Cannot load/empty Ps?',cl_id)), end
            % In the SW/SH we use a different set of offsets which are
            % independent of the spacecraft potential.
   			if caa_is_sh_interval
   				[dsiof_def, dam_def] = c_efw_dsi_off(t_int(1),cl_id,[]);
   			else
   				[dsiof_def, dam_def] = c_efw_dsi_off(t_int(1),cl_id,Ps);
   			end
   
   			[ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end
   %         [Ddsi, ok1] = caa_get(st, dt, cl_id, 'Ddsi?'); if ~ok1, Ddsi = dsiof_def; end
   			[ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end
   %			[Damp, ok2] = caa_get(st, dt, cl_id, 'Damp?'); if ~ok2, Damp = dam_def; end
   			
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
   		
%   		if length(Ddsi) == 1
%%   			dsi_str = sprintf('ISR2 offsets: dEx=%1.2f dEy=%1.2f dAmp=%1.2f',...
%            dsi_str = sprintf('Dsi offsets (ISR2): dEx=%1.2f dEy=%1.2f dAmp=%1.2f',...
%   				real(Ddsi(1)),imag(Ddsi(1)),Damp);
%   		else
%%   			dsi_str = 'ISR2 offsets';
%   			dsi_str = 'Dsi offsets';
%   			for in = 1:size(Ddsi,1)
%   				dsi_str = [dsi_str sprintf(' %s: dEx=%1.2f dEy=%1.2f,',...
%   					epoch2iso(Ddsi(in,1),1),real(Ddsi(in,2)),imag(Ddsi(in,2)))];
%   			end
%   			dsi_str = [dsi_str sprintf(' dAmp=%1.2f',Damp)];
%   		end
   		
   		if length(Ddsi) == 1
   		   dsi_str = sprintf('%s/%s ISR2 offsets: dEx=%1.2f dEy=%1.2f, dAmp=%1.2f',...
               epoch2iso(t_int(1),1),epoch2iso(t_int(2),1),real(Ddsi),imag(Ddsi),Damp);
         else
            for in = 1:(size(Ddsi, 1)-1)
               dsi_str = sprintf('%s/%s ISR2 offsets: dEx=%1.2f dEy=%1.2f, dAmp=%1.2f',...
                  epoch2iso(Ddsi(in,1),1),epoch2iso(Ddsi(in+1,1),1),real(Ddsi(in,2)),imag(Ddsi(in,2)),Damp);
               result_com{end+1} = dsi_str;
   	         irf_log('calb',dsi_str)
            end
            dsi_str = sprintf('%s/%s ISR2 offsets: dEx=%1.2f dEy=%1.2f, dAmp=%1.2f',...
               epoch2iso(Ddsi(end,1),1),epoch2iso(t_int(2),1),real(Ddsi(end,2)),imag(Ddsi(end,2)),Damp);
         end   
         result_com{end+1} = dsi_str;
   	   irf_log('calb',dsi_str)
   		clear Ddsi Damp
%   	   result_com{end+1} = dsi_str;
%   	   irf_log('calb',dsi_str)
   		
%   		dsc.com{dd} = ['Probe pair(s): ' E_info.probe ...
%   		   ' in subinterval: ' epoch2iso(t_int(1),1) '/' epoch2iso(t_int(2),1)];
%   		irf_log('calb',dsc.com)

%         result_com{end+1} = ['Probe pair(s): ' probe_info ...
%   		   ' in subinterval: ' epoch2iso(t_int(1),1) '/' epoch2iso(t_int(2),1)];
%   		irf_log('calb',result_com{end})
   		
   	end
   	
   	[ok,Del] = c_load('D?p12p34',cl_id); if ~ok, Del = [0 0]; end
   %   [Del, ok] = caa_get(st, dt, cl_id, 'D?p12p34'); if ~ok, Del = [0 0]; end
      if ~isreal(Del)
      	Del = imag(Del);
      	% offset is applied to p34
      	if strcmp(dsc.sen,'12') || strcmp(dsc.sen,'32'), Del = [0 0]; end
%      	del_str = sprintf('%s. p%s offset (ISR2): dEx=%1.2f dEy=%1.2f',...
%            dsi_str, '34', Del(1), Del(2));
      	del_str = sprintf('p%s offset (ISR2): dEx=%1.2f dEy=%1.2f', ...
      	   '34', Del(1), Del(2));
      else
      	% offset is applied to p12/32
      	if strcmp(dsc.sen,'34'), Del = [0 0]; end
%      	del_str = sprintf('%s. p%s offset (ISR2): dEx=%1.2f dEy=%1.2f',...
%            dsi_str, dsc.sen(1:2), Del(1), Del(2));
         del_str = sprintf('p%s offset (ISR2): dEx=%1.2f dEy=%1.2f',...
      		dsc.sen(1:2), Del(1), Del(2));
      end
      del_str = sprintf('%s/%s %s', ...
         epoch2iso(t_int(1),1), epoch2iso(t_int(2),1), del_str);
      result_com{end+1} = del_str;
      irf_log('calb',del_str)
   	
   	
   elseif lev==1 && ~isempty(regexp(caa_vs,'^P(12|32|34)?$','once'))
   	if ~isempty(data)
   		% convert mV/m back to V
   		if id==32, data(:,2) = data(:,2)*.0622;
           else data(:,2) = data(:,2)*.088;
   		end
   	end
   	
   % Combine ADC offsets from two probe pairs into one dataset:
   elseif strcmp(caa_vs, 'DER')
      if ~(isempty(data1) && isempty(data2))
   
         if length(probe_pairs) == 1   % Most likely p1 broken, so no 'Dadc?p12' product!
            single_pair_data = [data1 data2];
            start_time = min( single_pair_data(:,1) );
         else
            start_time = min( min(data1(:,1)), min(data2(:,1)) );
%            start_time = min( min( [data1(:,1) data2(:,1)] ) );
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
%            dsc.com = sprintf('Probe pairs used are p%i and p%i', probe_pairs);
            adc_str = sprintf('Probe pairs p%i,p%i', probe_pairs);
%            result_com{end+1} = [sprintf('Probe pairs: p%i,p%i', probe_pairs) ...
%   	   	   ' in subinterval: ' epoch2iso(t_int(1),1) '/' epoch2iso(t_int(2),1)];
         else
%            dsc.com = sprintf('Probe pair used is p%i', probe_pairs);
            adc_str = sprintf('Probe pair p%i', probe_pairs);
%            result_com{end+1} = [sprintf('Probe pair: p%i', probe_pairs) ...
%   	   	   ' in subinterval: ' epoch2iso(t_int(1),1) '/' epoch2iso(t_int(2),1)];
         end
         adc_str = sprintf('%s/%s %s', ...
            epoch2iso(t_int(1),1), epoch2iso(t_int(2),1), adc_str);
         result_com{end+1} = adc_str;
         irf_log('calb',adc_str)
         
         dsc.valtype = {dsc.valtype{:}, 'FLOAT'};
         dsc.sigdig = [dsc.sigdig 6];
         dsc.size = [dsc.size 1];
         
         clear start_time timestamp ind1 ind2 data_out
      
      else % isempty(data1) && isempty(data2)
         data = [];
      end
   end
   
   if isempty(result), result = data;
   else
	   if ~isempty(data)
		   t = result(:,1);
		   tapp = data(:,1);
		   
		   if tapp(1) <= t(end)
			   irf_log('proc',sprintf('Last point in data is %f seconds before first point in appended data',t(end)-tapp(1)))
			   irf_log('proc',sprintf('   Last point in data: %s',epoch2iso(t(end))))
			   irf_log('proc',sprintf('   Attempt to append interval %s to %s',epoch2iso(tapp(1)),epoch2iso(tapp(end))))
			   data = irf_tlim(data,tapp(1),t(end),1);
		   end
	   end
	   
       if ~isempty(data)
           result = [result; data];
       end
   end
       

   
%   if isempty(result), result = data;
%   else result = caa_append_data(result, data);    % NOTE: This will fill ALL columns with NaN unless
%                                                   % caa_fill_gaps is used BEFORE caa_identify_problems!
%   end

end   % for dd = 1:length(dirs)

data = result;
if exist('result_com') && length(result_com) > 0, dsc.com = result_com; end

cd(old_pwd)

% Check for non-monotonic time, and remove data within 2 HK packets (10.4s)
if (lev > 1) && ~isempty(data)
    indx=find(diff(data(:,1)) < 0.5e-3);
    if ~isempty(indx)
        irf_log('save',['WARNING: detected ' num2str(length(indx)) ' non-monotonic time stamp.'])
        for i=1:length(indx)
            irf_log('save',['Removing data near non-monotonic time stamp at ' epoch2iso(data(indx(i),1))])
            ii_left =find(data(:,1) < (data(indx(i)+1,1) - 10.4), 1, 'last' );
            ii_right=find(data(:,1) > (data(indx(i),1)   + 10.4), 1, 'first' );
            data(ii_left:ii_right,:)=NaN;             
        end
        indx=isfinite(data(:,1));
        data=data(indx,:);
    end
end

if isempty(data)
%   keyboard
   t_int_full = st + [0 dt];
   switch caa_vs
      case {'E', 'DER'}
         irf_log('save', sprintf('Saving empty interval %s/%s', ...
            epoch2iso(t_int_full(1),1), epoch2iso(t_int_full(2),1)) )
         v_size = 0;
         dsc.com = '';
         
%      case 'E'
%         dsc = c_desc(vs);
%         % Extend variable description to include the new columns bitmask and quality:
%         dsc.size = [dsc.size, 1, 1];
%         dsc.valtype = {dsc.valtype{:}, 'INT', 'INT'};
%         dsc.sigdig = [dsc.sigdig, 5, 1];
         
      otherwise
         status = 1; return
   end         
	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write to file
ext_s = '.cef';
% We have special names for CAA
DATASET_ID = irf_ssub(['C?_CP_EFW_L' num2str(lev) '_' caa_vs],cl_id);
if DELIVERY_TO_CAA
   file_name = ...
      [DATASET_ID '__' irf_fname(t_int_full,3) '_V' DATA_VERSION];
else
   file_name = ...
	   [DATASET_ID '__' irf_fname(t_int_full,2) '_V' DATA_VERSION];
end

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
['   ENTRY       =   ' epoch2iso(t_int_full(1)) '/' epoch2iso(t_int_full(2)) '\n']);
buf = sprintf('%s%s',buf,'END_META       =   FILE_TIME_SPAN\n');
buf = sprintf('%s%s',buf,'START_META     =   GENERATION_DATE\n');
buf = sprintf('%s%s',buf,'   VALUE_TYPE  =   ISO_TIME\n');
buf = sprintf('%s%s',buf,['   ENTRY       =   ' epoch2iso(date2epoch(nnow)) '\n']);
buf = sprintf('%s%s',buf,'END_META       =   GENERATION_DATE\n');
if strcmp(caa_vs, 'SFIT')
    if nanfill
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'P34 data only.' dsc.com ]);
    else
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'P' num2str(pnosfit) ' & P34 data.' dsc.com ]);
    end
else
    buf = pmeta(buf, 'FILE_CAVEATS', dsc.com);
end    

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
else
   [fid,msg] = fopen([file_name ext_s],'a');
   if fid < 0
	   irf_log('save',['problem opening CEF file: ' msg])
	   status = 1;
	   return
   end
   
   sta = fprintf(fid,'END_OF_DATA\n');
   fclose(fid);
   if sta<=0, irf_log('save','problem writing CEF tail'), status = 1; return, end 
   
   [sta, res] = unix(['gzip ' file_name ext_s]);
   if sta>0
      irf_log('save', ['Error gzipping output file: ' res])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dir_list = get_subdirs(iso_t, dt, cl_id)
% Find all subinterval dirs for spacecraft cl_id
% NOTE: For the export of 24-h files, this has been replaced
%       by caa_get_subdirs.m

SPLIT_INT = 3; % 3 hour subintervals
BASE_DIR = '/data/caa/l1';

[st,dt] = irf_stdt(iso_t,dt);

t = fromepoch(st);
t0 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);
t = fromepoch(st+dt);
t1 = toepoch([t(1) t(2) t(3) fix(t(4)/SPLIT_INT)*SPLIT_INT 0 0]);
if t1>=st+dt, t1 = t1 - SPLIT_INT*3600; end

old_pwd = pwd;
mode_list = [];
for t=t0:SPLIT_INT*3600:t1
	y = fromepoch(t);
	main_int = [BASE_DIR '/' num2str(y(1)) '/' irf_fname(t) '/C' num2str(cl_id)];
	if ~exist(main_int,'dir'), continue, end
	
	cd(main_int)
	d = dir('*_*');
	if isempty(d), continue, end
	good_dir = {};
	for j=1:length(d)
		if ~d(j).isdir, continue, end
		if caa_is_valid_dirname(d(j).name), good_dir = {good_dir{:} d(j).name}; end
	end
	if isempty(good_dir), continue, end
end
dir_list = good_dir;
	
%	for j=1:length(good_dir)
%		subdir = [main_int '/' good_dir{j}];
%		cd(subdir)
%		[st_t,dt_tmp] = caa_read_interval();
%		if isempty(st_t), continue, end