function status = caa_export_cef(lev,caa_vs,cl_id,QUALITY,DATA_VERSION,sp,st,dt)
%CAA_EXPORT_CEF  export data to CAA CEF files
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

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

status = 0;

% This must be changed when we do any major changes to our processing software
% Version 5: new algorithm to compute ScPot from several probes, median
EFW_DATASET_VERSION = '5';

if nargin~=8 && nargin~=3, error('3 or 8 input parameters are needed'); end
if nargin ==3
    QUALITY = 3;
    DATA_VERSION = '00';
    DELIVERY_TO_CAA = 0;    % testing mode
else
    DELIVERY_TO_CAA = 1;    % Changes file name format to new CAA daily-file-format
end
if cl_id<=0 || cl_id>4, error('CL_ID must be 1..4'), end
if lev<1 || lev>3, error('LEV must be 1,2 or 3'), end

DATASET_DESCRIPTION_PREFIX = '';
EOR_MARKER = '$';
FILL_VAL = -1.0E9;
PROCESSING_LEVEL='Calibrated';


old_pwd = pwd;

if 1
if DELIVERY_TO_CAA   % Check that files are midnight-to-midnight
   st_temp = fromepoch(st);
   if dt ~= 24*3600 || st_temp(4) ~= 00
      error('Incorrect format for daily files. Check parameters!')
   end
end
end
ibsave=false;
badibfound=false;
result = [];
result_com = {};


if lev==1
	if strcmp(caa_vs, 'IB')
        vs = irf_ssub('IB?',cl_id);
		v_size = 8;
        alldi='';
    elseif regexp(caa_vs,'^P(1|2|3|4|12|32|34)?$')
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
		if lev==2 || lev==3
            % Fake for c_desc only. No data variable in .mat files
			vs = irf_ssub('P?',cl_id);
			v_size = 5;
		else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'E'
        if lev==2
            vs = irf_ssub('diE?p1234',cl_id);
            v_size = 3;    % Previously 1. (ML)
        elseif lev==3
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
	    vs_DER = sprintf('Dadc%.0fp?', cl_id);
	    probe_pair_list = [12, 32, 34];
		v_size = 2;
	case 'SFIT'
		if lev==3
            % Fake for c_desc only. No data variable in .mat files
			vs = irf_ssub('SFIT?',cl_id);
			v_size = 4;
            nanfill = -1;
            pnosfit = 0;
        else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'PB'
		if lev==2
            % Fake for c_desc only. No data variable in .mat files
			vs = irf_ssub('PB?',cl_id);
			v_size = 5;
        else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'BB'
		if lev==2
            % Fake for c_desc only. No data variable in .mat files
			vs = irf_ssub('BB?',cl_id);
			v_size = 5;
        else
			disp('not implemented'), cd(old_pwd), return
		end
	case 'EB'
		if lev==2
            % Fake for c_desc only. No data variable in .mat files
			vs = irf_ssub('EB?',cl_id);
			v_size = 4;
        else
			disp('not implemented'), cd(old_pwd), return
		end
	otherwise
		error(['Unknown variable ' caa_vs])
	end
end

if DELIVERY_TO_CAA, dirs = caa_get_subdirs(st, dt, cl_id);
else
    dirs = {'.'};
    [iso_st,dt] = caa_read_interval();
    st = iso2epoch(iso_st);
end
for dd = 1:length(dirs)
   d = dirs{dd};
   cd(d);
   
   % Set up spin fit related information.
   % Probe pair used can vary for each subinterval!
   if strcmp(caa_vs, 'E')
      [sfit_probe,flag_lx,probeS] = caa_sfit_probe(cl_id);
      irf_log('proc',sprintf('using %s',probeS))
      if lev == 3
          if sfit_probe == 320
              sfit_probe = 32;
          end
          if sfit_probe == 340
              sfit_probe = 34;
          end
          if sfit_probe == 420
              sfit_probe = 42;
          end
          if flag_lx, vs = irf_ssub('diELXs?p!',cl_id,sfit_probe);
          else vs = irf_ssub('diEs?p!',cl_id,sfit_probe);
          end
      end
   end

   % Load data
   if strcmp(caa_vs, 'DER')
      vs = vs_DER;
      [ok, data] = c_load(vs,cl_id,'res',probe_pair_list);   % Try loading data for all probe pairs.

      probe_pairs = probe_pair_list(logical(ok));   % Keep list of probe pairs actually loaded.
      if numel(probe_pairs) > 0
         vs = irf_ssub(vs, probe_pairs(1));
      end
   elseif strcmp(caa_vs, 'P')      
      if ~exist('c_ct','var')
         global c_ct % includes aspoc active values
      end
      if isempty(c_ct)
         c_ctl('load_aspoc_active', [c_ctl('get', 5, 'data_path') '/caa-control']);
         global c_ct
      end
      ASPOC = c_ct{1,cl_id}.aspoc;

      if lev == 2
          pvar = 'P?';
      elseif lev == 3
          pvar = 'Ps?';
      end
      
      [ok,probe_info,msg] = c_load([ pvar '_info'],cl_id);
      if ~ok || isempty(probe_info)
         irf_log('load',msg)
         data = [];
      else
         [ok,data,msg] = c_load(pvar,cl_id);
         if ~ok || isempty(data)
            irf_log('load',msg)
            data = [];
         else
            % Extend data array to accept probe#, aspoc, bitmask and quality (4 new columns at the end)
            if size(data,1) == 1    % Fix short data
                data=[data;[data(1,1)+4 NaN]];
                irf_log('proc','short data padded');
            end
            dsize=size(data, 1);
            data = [data zeros(dsize, 4)]; % add columns: probe# aspoc bitmask quality
            data(:, 3) = probe_info.probe; % Set probe#.

            if isempty(ASPOC)
                irf_log('proc','no ASPOC active data');
            else
                for i=1:size(ASPOC,1)
                    %i
                    if data(1,1)>ASPOC(i,1) && data(1,1)>ASPOC(i,2) %  too early: next
                        continue;
                    end
                    if data(dsize,1)<ASPOC(i,1) && data(dsize,1)<ASPOC(i,2) % too late: stop
                        break;
                    end
                    irf_log('proc','marking ASPOC active');
                    for j=1:dsize
                       if data(j,1)>=ASPOC(i,1) && data(j,1)<=ASPOC(i,2)
                            data(j,4)=1;
                        end
                    end
                end
            end
            data(:, end) = QUALITY;        % Default quality column to best quality, i.e. good data/no problems.
            quality_column = size(data, 2);
            bitmask_column = quality_column - 1;

            % Identify and flag problem areas in data with bitmask and quality factor:
            data = caa_identify_problems(data, lev, num2str(probe_info.probe), cl_id, bitmask_column, quality_column, 1);

         end
      end
   elseif strcmp(caa_vs, 'SFIT') 
          no_p12 = 0; no_p34 = 0;
          % Load P34
          [ok,spf34,msg] = c_load('diEs?p34',cl_id);
          if ~ok || isempty(spf34)
              no_p34 = 1;
              irf_log('load',msg)
          end
          % Load P12/32
          if exist('./mEDSI.mat')
              pnosfit = 12;
              ret=whos('-file','./mEDSI.mat',irf_ssub('diEs?p!',cl_id,pnosfit));
              if isempty(ret)
                 pnosfit = 32;
              end
              [ok,spfD,msg] = c_load(irf_ssub('diEs?p!',cl_id,pnosfit));
              if ~ok || isempty(spfD)
                 no_p12 = 1; % No P12/32 data
                 irf_log('load',msg)    
              end
          else
              no_p12 = 1;
          end
          if no_p12 && no_p34
             data = [];
          elseif no_p12
             % save time, NaN(fillval) and p34 spin-fit (B C sdev)
             nanfill = 0;
             data=[spf34(:,1) NaN(size(spf34,1),3) spf34(:,2:3) spf34(:,5)];
          elseif no_p34
             nanfill = 1;
             % save time, p12/32 spin-fit (B C sdev) and NaN(fillval)  
             data=[spfD(:,1) spfD(:,2:3) spfD(:,5) NaN(size(spfD,1),3)];
          else
             % save time, p12/32 and p34 spin-fit: B C sdev
             [ok,Del,msg] = c_load('D?p12p34',cl_id);
             if ~ok || isempty(Del)
                 irf_log('load',msg)
                 data = []; continue
             end;
             if imag(Del(1)) ~= 0 || imag(Del(2)) ~= 0
                 irf_log('load','Info: Imaginary delta offset.');
             end
             spfD(:,2:3)=spfD(:,2:3)+ones(size(spfD,1),1)*Del;
             s34=size(spf34(:,1),1);
             sd=size(spfD(:,1),1);
             if s34 > sd
                found=0;
                for j=1:sd
                    for i=1:s34
                        if spf34(i,1) == spfD(j,1)
                            found=i;
                            pos=j;
                            break;
                        end
                    end
                    if found
                        break;
                    end
                end
                if found
                   if pos > 1
                      spf34=[NaN(pos-1,5);spf34];
                      irf_log('proc','Info: pos>1 sd');
                   end
                   spfD=[NaN(found-1,5);spfD];
                   if s34-found-sd >= 0
                       spfD=[spfD;NaN(s34-found-sd+pos,5)];
                   end
                else
                   if spf34(1,1) < spfD(1,1)
                       spf34=[spf34;NaN(sd,5)];
                       spfD=[NaN(s34,5);spfD];
                   else
                       spf34=[NaN(sd,5);spf34];
                       spfD=[spfD;NaN(s34,5)];
                   end    
                end
                data=[spf34(:,1) spfD(:,2:3) spfD(:,5) spf34(:,2:3) spf34(:,5)];
             elseif s34 < sd
                found=0;
                for j=1:s34
                    for i=1:sd
                        if spf34(j,1) == spfD(i,1)
                            found=i;
                            pos=j;
                            break;
                        end
                    end
                    if found
                        break;
                    end
                end
                if found
                   if pos > 1
                      spfD=[NaN(pos-1,5);spfD];
                      irf_log('proc','Info: pos>1 s34');
                   end
                   spf34=[NaN(found-1,5);spf34];
                   if sd-found-s34 >= 0
                       spf34=[spf34;NaN(sd-found-s34+pos,5)];
                   end
                else
                   if spf34(1,1) < spfD(1,1)
                       spf34=[spf34;NaN(sd,5)];
                       spfD=[NaN(s34,5);spfD];
                   else
                       spf34=[NaN(sd,5);spf34];
                       spfD=[spfD;NaN(s34,5)];
                   end    
                end
                data=[spfD(:,1) spfD(:,2:3) spfD(:,5) spf34(:,2:3) spf34(:,5)];
             else
                data=[spf34(:,1) spfD(:,2:3) spfD(:,5) spf34(:,2:3) spf34(:,5)];
             end
          end
   elseif strcmp(caa_vs, 'IB')
       ok=0;
       if lev==1
         mfn='./mEFWburstTM.mat'; % For tm
         if exist(mfn,'file')
           r=load(mfn);
           di=eval(irf_ssub('r.ib?_info',cl_id));
           if isempty(alldi)
                alldi=di;
           else
                alldi=[alldi ', ' di];
           end
           data=eval(irf_ssub('r.iburst?',cl_id));
           ok=1;
         else
           data=[];
         end
       end
   elseif strcmp(caa_vs, 'PB')
       if lev==2
         if ~exist('c_ct','var')
            global c_ct % includes aspoc active values
         end
         if isempty(c_ct)
            c_ctl('load_aspoc_active', [c_ctl('get', 5, 'data_path') '/caa-control']);
%            c_ctl('load_aspoc_active');
            global c_ct
         end
         ASPOC = c_ct{1,cl_id}.aspoc;

         % Export empty file if TM, BB or EB data exist
        % Check for TM data
         mfn='./mEFWburstTM.mat'; % For tm
         if exist(mfn,'file')
             r=load(mfn);
             data=eval(irf_ssub('r.iburst?',cl_id));
             if ~isempty(data)
                 ibsave=true;
                 clear data;
             end
         else
              % Check for B data
             mfn='./mEFWburst.mat';
             if exist(mfn,'file')
                 bsc=load(mfn);
                 finbsc=fieldnames(bsc);
                 fnl=size(finbsc,1);
                 found=0;
                 for bscix=1:fnl % find BSC data
        %               finbsc{bscix}
                   if length(finbsc{bscix})<5
                       continue;
                   end
                   if strcmp(finbsc{bscix}(1:5),'diBSC')
                     ibsave=true;
                     break;
                   end
                 end
             end
             if ~ibsave
                 % Check for E data
                 mfn='./mEFWburst.mat'; % For E
                 if exist(mfn,'file')
                     e=load(mfn);
                     fine=fieldnames(e);
                     fnl=size(fine,1);
                     found=0;
                     for eix=1:fnl % find dibE data
        %               fine{eix}
                       if length(fine{eix})<4
                           continue;
                       end
                       if strcmp(fine{eix}(1:4),'dibE') && ~strcmp(fine{eix}(end-3:end),'info')
                         ibsave=true;
                         break;
                       end
                     end
                 end
             end
         end

         pvar='bP?';
         [ok,probe_info,msg] = c_load([ pvar '_info'],cl_id);
         if ~ok || isempty(probe_info) % Check for no IB data
%            irf_log('load',msg)
            data = [];
         else
            [ok,data,msg] = c_load(pvar,cl_id);
            if ~ok || isempty(data)
                irf_log('load',msg)
                data = [];
            else
                % Extend data array to accept probe#, aspoc, bitmask and quality (4 new columns at the end)
                if size(data,1) == 1    % Fix short data
                    data=[data;[data(1,1)+4 NaN]];
                    irf_log('proc','short data padded');
                end
                dsize=size(data, 1);
                data = [data zeros(dsize, 4)]; % add columns: probe# aspoc bitmask quality
                data(:, 3) = probe_info.probe; % Set probe#.

                if isempty(ASPOC)
                    irf_log('proc','no ASPOC active data');
                else
                    for i=1:size(ASPOC,1)
                        %i
                        if data(1,1)>ASPOC(i,1) && data(1,1)>ASPOC(i,2) %  too early: next
                            continue;
                        end
                        if data(dsize,1)<ASPOC(i,1) && data(dsize,1)<ASPOC(i,2) % too late: stop
                            break;
                        end
                        irf_log('proc','marking ASPOC active');
                        for j=1:dsize
                           if data(j,1)>=ASPOC(i,1) && data(j,1)<=ASPOC(i,2)
                                data(j,4)=1;
                            end
                        end
                    end
                end
                data(:, end) = QUALITY;        % Default quality column to best quality, i.e. good data/no problems.
                quality_column = size(data, 2);
                bitmask_column = quality_column - 1;

                % Identify and flag problem areas in data with bitmask and quality factor:
                data = caa_identify_problems(data, lev, num2str(probe_info.probe), cl_id, bitmask_column, quality_column, 2);
            end
         end
       end
   elseif strcmp(caa_vs, 'BB')
       if lev==2
         if ~exist('c_ct','var')
            global c_ct % includes aspoc active values
         end
         if isempty(c_ct)
            c_ctl('load_bad_ib', [c_ctl('get', 5, 'data_path') '/caa-control'])
%            c_ctl('load_bad_ib');
            global c_ct
         end

         % Export empty file if TM, PB or EB data exist. Ex 110227205958we.03
         % Check for TM data
         mfn='./mEFWburstTM.mat'; % For tm
         if exist(mfn,'file')
             r=load(mfn);
             data=eval(irf_ssub('r.iburst?',cl_id));
             if ~isempty(data)
                 ibsave=true;
             end
         else
             % Check for P data
             pvar='bP?';
             [ok,data,msg] = c_load(pvar,cl_id);
             if ok && ~isempty(data)
                ibsave=true;
             else
                 % Check for E data
                 mfn='./mEFWburst.mat';
                 if exist(mfn,'file')
                     e=load(mfn);
                     fine=fieldnames(e);
                     fnl=size(fine,1);
                     found=0;
                     for eix=1:fnl % find dibE data
        %               fine{eix}
                       if length(fine{eix})<4
                           continue;
                       end
                       if strcmp(fine{eix}(1:4),'dibE') && ~strcmp(fine{eix}(end-3:end),'info')
                         ibsave=true;
                         break;
                       end
                     end
                 end
             end
         end
         data=[];
         ok=0;
         mfn='./mEFWburst.mat'; % Read Bx By Bz
         if exist(mfn,'file')
             bsc=load(mfn);
             finbsc=fieldnames(bsc);
             fnl=size(finbsc,1);
             found=0;
             for bscix=1:fnl % find BSC data
%               finbsc{bscix}
               if length(finbsc{bscix})<5
                   continue;
               end
               if strcmp(finbsc{bscix}(1:5),'diBSC')
                 found=1;
                 break;
               end
             end
             if found
               irf_log('proc','BSC burst data found');
               data=eval(['bsc.' finbsc{bscix}]);
               ok=1;
               % Extend data array to accept probe#, aspoc, bitmask and quality (4 new columns at the end)
               if size(data,1) == 1    % Fix short data
                   data=[data;[data(1,1)+4 NaN NaN NaN]];
                   irf_log('proc','short data padded');
               end
               dsize=size(data, 1);
               data = [data zeros(dsize, 2)]; % add columns: bitmask quality

               data(:, end) = QUALITY;        % Default quality column to best quality, i.e. good data/no problems.
               quality_column = size(data, 2);
               bitmask_column = quality_column - 1;

               % Identify and flag problem areas in data with bitmask and quality factor:
               data = caa_identify_problems(data, lev, num2str(1), cl_id, bitmask_column, quality_column, 4);
             else
                irf_log('load','No diBSC data matrix in mEFWburst.mat. No BB.');
             end
         end
       end
   elseif strcmp(caa_vs, 'EB')
       if lev==2
         if ~exist('c_ct','var')
            global c_ct % includes aspoc active values
         end
         if isempty(c_ct)
            c_ctl('load_bad_ib', [c_ctl('get', 5, 'data_path') '/caa-control'])
%            c_ctl('load_bad_ib');
            global c_ct
         end

         % Export empty file if TM, PB or BB data exist. Ex 110227205958we.03
         % Check for TM data
         mfn='./mEFWburstTM.mat'; % For tm
         if exist(mfn,'file')
             r=load(mfn);
             data=eval(irf_ssub('r.iburst?',cl_id));
             if ~isempty(data)
                 ibsave=true;
             end
         else
             % Check for P data
             pvar='bP?';
             [ok,data,msg] = c_load(pvar,cl_id);
             if ok && ~isempty(data)
                ibsave=true;
             else
                 % Check for B data
                 mfn='./mEFWburst.mat';
                 if exist(mfn,'file')
                     bsc=load(mfn);
                     finbsc=fieldnames(bsc);
                     fnl=size(finbsc,1);
                     found=0;
                     for bscix=1:fnl % find BSC data
        %               finbsc{bscix}
                       if length(finbsc{bscix})<5
                           continue;
                       end
                       if strcmp(finbsc{bscix}(1:5),'diBSC')
                         ibsave=true;
                         break;
                       end
                     end
                 end
             end
         end
         data=[];
         ok=0;
         mfn='./mEFWburst.mat'; % For E
         if exist(mfn,'file')
             e=load(mfn);
             fine=fieldnames(e);
             fnl=size(fine,1);
             found=0;
             for eix=1:fnl % find dibE data
%               fine{eix}
               if length(fine{eix})<4
                   continue;
               end
               if strcmp(fine{eix}(1:4),'dibE') && ~strcmp(fine{eix}(end-3:end),'info')
                 found=1;
                 break;
               end
             end
             if found
               irf_log('proc','E burst data found');
               [ok,probe_info,msg] = c_load([ fine{eix} '_info']);
               data=eval(['e.' fine{eix}]);
               ok=1;
               % Extend data array to accept bitmask and quality (1 new column at the end 1 reused)
               if size(data,1) == 1    % Fix short data
                   data=[data;[data(1,1)+4 NaN NaN NaN]];
                   irf_log('proc','short data padded');
               end
               dsize=size(data, 1);
               data = [data zeros(dsize, 1)]; % add columns: bitmask quality
%               data(:, 4) = str2num(probe_info.probe); % Set probe#.

               data(:, end) = QUALITY;        % Default quality column to best quality, i.e. good data/no problems.
               quality_column = size(data, 2);
               bitmask_column = quality_column - 1;

               % Identify and flag problem areas in data with bitmask and quality factor:
               data = caa_identify_problems(data, lev, probe_info.probe, cl_id, bitmask_column, quality_column, 3);
             else
                irf_log('load','No dibE data matrix in mEFWburst.mat. No EB.');
             end
         end
       end
   else
      [ok,data] = c_load(vs);
   end
   if (all(~ok) || isempty(data))
        if ~regexp(caa_vs,'^(I|P|E|B)B$') 
            irf_log('load', ['No ' vs]);
        end
        cd(old_pwd)
        continue
   end
   d_info = []; ok = 0;
   try
      if ~strcmp(caa_vs, 'SFIT') && ~regexp(caa_vs,'^(I|P|E|B)B$') 
        [ok, d_info] = c_load([vs '_info'],'var');
      else
   	    d_info = []; ok = 0;
      end
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
       if regexp(caa_vs,'^(I|P|E|B)B$') % iburst can be only partly inside time frame
           t_int(1)=t_int(1)-120;
           t_int(2)=t_int(2)+300;
       end
           
         % Limit data to both time interval given as input, and time interval read from file.
   	   data = irf_tlim(data,t_int_full);  % NOTE: Superfluous when using caa_get above. (ML)
   	   data = irf_tlim(data, t_int);
   	   if isempty(data)
   		   irf_log('save', 'Saving empty subinterval')
       end
       % Check for bad iburst file
%if 0
        if regexp(caa_vs,'^(P|E|B)B$')   % working on multiple iburst in 24h
           if ~exist('c_ct','var')
                global c_ct % includes bad ib files
           end
           if isempty(c_ct{1,1}.badib)
                c_ctl('load_bad_ib', [c_ctl('get', 5, 'data_path') '/caa-control']);
                global c_ct
           end
           badiburst = c_ct{1,cl_id}.badib;
           tint = irf_tlim(badiburst,t_int_full);
           if ~isempty(tint) % remove iburst
               for i=1:size(tint,1)
                    if ~isempty(data)
                        sizemem=size(data,1);
                        data = irf_tlim(data,tint(i)-90,tint(i)+240,1); % remove -1.5 min to +4 min
                        if size(data,1)~=sizemem
                            badibfound=true;
                            irf_log('save', [caa_vs num2str(cl_id) ' L2 iburst removed. Data marked as bad.'])
                        end

                    end
               end
           end
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
if ~isempty(data)
    tdiff=0.5e-3;
    if regexp(caa_vs,'^(I|P|E|B)B$') 
        tdiff=1e-5;
    end
    indx=find(diff(data(:,1)) < tdiff);
    if ~isempty(indx)
        irf_log('save',['WARNING: detected ' num2str(length(indx)) ' non-monotonic time stamp.'])
        for i=1:length(indx)
            irf_log('save',['Removing data near non-monotonic time stamp at ' epoch2iso(data(indx(i),1))])
            ii_left =find(data(:,1) < (data(indx(i)+1,1) - 10.4), 1, 'last' );
			if isempty(ii_left), ii_left=1; end
            ii_right=find(data(:,1) > (data(indx(i),1)   + 10.4), 1, 'first' );
			if isempty(ii_right), ii_right=length(data(:,1)); end
            data(ii_left:ii_right,:)=NaN;             
        end
        indx=isfinite(data(:,1));
        data=data(indx,:);
        if strcmp(caa_vs,'IB') && isempty(data)
            ibsave=true;
        end 
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
         dsc.com = '';
%         status = 1; return
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
    if pnosfit == 0
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'No data.' dsc.com ]);
    elseif nanfill == 0
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'P34 data only.' dsc.com ]);
    elseif nanfill == 1
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'P' num2str(pnosfit) ' data only.' dsc.com ]);
    else
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'P' num2str(pnosfit) ' & P34 data.' dsc.com ]);
    end
elseif strcmp(caa_vs, 'IB')
    if lev==1
        if isempty(data)
            alldi='na';
        end
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'Data order: ' alldi ' ' dsc.com ]);
    end
elseif ~isempty(regexp(caa_vs,'^(P|E|B)B$')) && ibsave && isempty(data)
        buf = pmeta(buf, 'FILE_CAVEATS', [ 'No iburst ' caa_vs ' data. ' dsc.com ]);    
else
    buf = pmeta(buf, 'FILE_CAVEATS', dsc.com);
end    

buf = sprintf('%s%s',buf,'!\n');
buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,'!                       Data                          !\n');
buf = sprintf('%s%s',buf,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
buf = sprintf('%s%s',buf,'DATA_UNTIL = "END_OF_DATA"\n');

if (~strcmp(caa_vs, 'IB') && ~strcmp(caa_vs, 'PB') && ~strcmp(caa_vs, 'EB') && ~strcmp(caa_vs, 'BB')) || ~isempty(data) || ibsave || badibfound
    [fid,msg] = fopen([file_name ext_s],'w');
    if fid < 0
        irf_log('save',['problem opening CEF file: ' msg])
        status = 1;
        return
    end

    sta = fprintf(fid,buf);
    fclose(fid);
    if sta<=0, irf_log('save','problem writing CEF header'), status = 1; return, end
end
if ~isempty(data)
	n_col = size(data,2) -1; % number of data columns - time
	for j=1:n_col
		ii = find(abs(data(:,j+1)) > 1e8);
		if ~isempty(ii)
			irf_log('save',['WARNING: detected ' num2str(length(ii)) ' wildly out of range data points.'])
			data(ii,j+1) = FILL_VAL;
		end
		ii = find(isnan(data(:,j+1))) ;
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
elseif ~isempty(regexp(caa_vs,'^(I|P|E|B)B$')) && ~ibsave && ~badibfound 
   irf_log('proc','Will not export empty internal burst IB, PB, EB or BB files')
else
   disp(['Filename : ' file_name ext_s ' (Empty)' ]);
   [fid,msg] = fopen([file_name ext_s],'a');
   if fid < 0
	   irf_log('save',['problem opening CEF file: ' msg])
	   status = 1;
	   return
   end
   
   sta = fprintf(fid,'END_OF_DATA\n');
   fclose(fid);
   if sta<=0, irf_log('save','problem writing CEF tail'), status = 1; return, end 

   if exist([file_name ext_s '.gz'],'file')
       delete([file_name ext_s '.gz']);
   end
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

