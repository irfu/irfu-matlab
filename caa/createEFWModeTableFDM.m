function [t_start_save,t_stop_save,fdm_save]=createEFWModeTableFDM(db_s,start_epoch,delta_t,cli,var_list_s)
%createEFWModeTableFDM reads EFW FDM and extracts diffs
% [t_start_save,t_stop_save,fdm_save]=
% 	createEFWModeTableFDM(db_s,start_epoch,dt,cl_id,var_list_s)
%
% db_s	ISDAT database string ('disco:10')
% var_list_s is in form 'px|r|ss'
% where:
% px	Burst playback
% bbb	Burst internal state
% r		Whisper pulses present
% eeee	E/D mode
% ss	Sampling mode
% tm	Tape mode
%
% $Revision$  $Date$
%

DEBUG=0;

dt_max = 60; % maximum interval is 1 min

stop_epoch = start_epoch + delta_t;

if DEBUG
	s = fromepoch(start_epoch);
	e = fromepoch(stop_epoch);
	disp(sprintf('Processing EFW modes from %d-%d-%d %d:%d:%2.3f to %d-%d-%d %d:%d:%2.3f', ...
	s(1),s(2),s(3),s(4),s(5),s(6),e(1),e(2),e(3),e(4),e(5),e(6)))
end

dt = delta_t;
if dt > dt_max, dt = dt_max; end
t_cur = start_epoch;

% storage variables
t_start_save = [];
t_stop_save = [];
fdm_save = {};
n_good = 0;
n_jumpy = 0;
l = 0;

% temporal variables
t_fdm_last = [];
fdm_last = [];

prev_ok = 0; % indicate if last time interval contained data
no_data = 1; % indicate if last time interval contained no data
t_end_last_interv = []; % end of last time interval;

var_list = tokenize(var_list_s,'|');

while t_cur < stop_epoch
	[t_fdm,fdm_is] = ISGet(db_s,t_cur,dt,cli,'efw','FDM');
	
	t_cur = t_cur + dt;
	if DEBUG
		disp(sprintf('%2.1f %% done', ...
		100*(t_cur - start_epoch)/(stop_epoch - start_epoch)))
	end

	if ~isempty(t_fdm);
		% this info is needed if next interval has no data
		prev_ok = 1;
		t_end_last_interv = t_fdm(end);
	
		% convert to bits
		s_pbbbxrwc = dec2bin(fdm_is(1,:),8);
		s_ssqiiiii = dec2bin(fdm_is(2,:),8);
		s_mmmmeeee = dec2bin(fdm_is(4,:),8);
		fdm_tm=fdm_is(5,:)';
		
		clear fdm_is
		
		fdm_px = str2num( s_pbbbxrwc(:,1) ); 
		fdm_bbb = str2num( s_pbbbxrwc(:,2) ); 
		fdm_bbb(:,2) = str2num( s_pbbbxrwc(:,3) ); 
		fdm_bbb(:,3) = str2num( s_pbbbxrwc(:,4) ); 
		fdm_px(:,2) = str2num( s_pbbbxrwc(:,5) );
		fdm_r = str2num( s_pbbbxrwc(:,6) );
		fdm_s = str2num( s_ssqiiiii(:,1) );
		fdm_s(:,2) = str2num( s_ssqiiiii(:,2) );
		fdm_e = str2num( s_mmmmeeee(:,5) );
		fdm_e(:,2) = str2num( s_mmmmeeee(:,6) );
		fdm_e(:,3) = str2num( s_mmmmeeee(:,7) );
		fdm_e(:,4) = str2num( s_mmmmeeee(:,8) );
		
		clear s_pbbbxrwc s_ssqiiiii s_mmmmeeee
		
		if no_data % if we just started or last interval had no data
			t_fdm_last = t_fdm(1);
			for j=1:length(var_list)
				if strcmp(var_list(j),'px')
					fdm_last_px = fdm_px(1,:);
				elseif strcmp(var_list(j),'bbb')
					fdm_last_bbb = fdm_bbb(1,:);
				elseif strcmp(var_list(j),'r')
					fdm_last_r = fdm_r(1,:);
				elseif strcmp(var_list(j),'s')
					fdm_last_s = fdm_s(1,:);
				elseif strcmp(var_list(j),'e')
					fdm_last_e = fdm_e(1,:);
				elseif strcmp(var_list(j),'tm')
					fdm_last_tm = fdm_tm(1)
				else
					error('unknown variable')
				end
			end
			no_data = 0;
		end
		
		for j=1:length(var_list)
			eval(['ii{j} = find_diff(fdm_' var_list{j} ',fdm_last_' var_list{j} ');']);
		end
		
		if isempty(ii)
			continue % no changes (should be a typical situation)
		end
		
		ii_all = [];
		for j=1:length(var_list)
			if ~isempty(ii{j}) 
				n = length(ii{j});
				ii_all(end+1:end+n) = ii{j};
			end
		end
		ii_all = rmDoubleIndex(ii_all);
		n_diff = length(ii_all);
		
		for jj = 1:n_diff
			l = length(t_start_save);
			
			% save a good interval
			t_start_save(l+1) = t_fdm_last;
			% if the change happened exactly at the boundary between intervals
			if ii_all(jj)-1 >= 1, t_stop_save(l+1) = t_fdm(ii_all(jj)-1);
			else t_stop_save(l+1) = t_end_last_interv; 
			end

			for j=1:length(var_list)
				eval(['if length(fdm_last_' var_list{j} ')>1,fdm_save_' var_list{j} '(l+1,:)=fdm_last_' var_list{j} '; else, fdm_save_' var_list{j} '(l+1)=fdm_last_' var_list{j} '; end ']);
			end
			%fdm_save{l+1} = fdm_last;
			n_good = n_good + 1;
			
			t_fdm_last = t_fdm(ii_all(jj));
			for j=1:length(var_list)
				eval(['fdm_last_' var_list{j} '= fdm_' var_list{j} '(ii_all(jj),:);']); 
			end
				
			if DEBUG	
				s = fromepoch(t_start_save(l+1));
				e = fromepoch(t_stop_save(l+1));
				disp(sprintf('Saving good from %d-%d-%d %d:%d:%2.3f to %d-%d-%d %d:%d:%2.3f', ...
				s(1),s(2),s(3),s(4),s(5),s(6),e(1),e(2),e(3),e(4),e(5),e(6)))
			end
		end
	
	else % no data

		no_data = 1;
		
		s = fromepoch(t_cur - dt);
		disp(sprintf('No data at %d-%d-%d %d:%d:%2.3f',...
		s(1),s(2),s(3),s(4),s(5),s(6)))
		
		if prev_ok
			l = length(t_start_save);
			
			% save the good interval
			t_start_save(l+1) = t_fdm_last;
			t_stop_save(l+1) = t_end_last_interv;

			%fdm_save{l+1} = fdm_last;
			for j=1:length(var_list)
				eval(['if length(fdm_last_' var_list{j} ')>1,fdm_save_' var_list{j} '(l+1,:)=fdm_last_' var_list{j} '; else, fdm_save_' var_list{j} '(l+1)=fdm_last_' var_list{j} '; end ']);
			end
			n_good = n_good + 1;
		
			if DEBUG	
				s = fromepoch(t_start_save(l+1));
				e = fromepoch(t_stop_save(l+1));
				disp(sprintf('Saving good from %d-%d-%d %d:%d:%2.3f to %d-%d-%d %d:%d:%2.3f', ...
				s(1),s(2),s(3),s(4),s(5),s(6),e(1),e(2),e(3),e(4),e(5),e(6)))
			end
			
			%fdm_last = [];
			for j=1:length(var_list)
				eval(['fdm_last_' var_list{j} '=[];']);
			end
			prev_ok = 0;
		end
	end
end

if prev_ok %last interval was also good
if (t_fdm_last < t_fdm(end))
	l = length(t_start_save);
	
	% save the good interval		
	t_start_save(l+1) = t_fdm_last;
	t_stop_save(l+1) = t_fdm(end);

	%fdm_save{l+1} = fdm_last;
	for j=1:length(var_list)
		eval(['if length(fdm_last_' var_list{j} ')>1,fdm_save_' var_list{j} '(l+1,:)=fdm_last_' var_list{j} '; else, fdm_save_' var_list{j} '(l+1)=fdm_last_' var_list{j} '; end ']);
	end
	n_good = n_good + 1;
					
	if DEBUG	
		s = fromepoch(t_start_save(l+1));
		e = fromepoch(t_stop_save(l+1));
		disp(sprintf('Saving good from %d-%d-%d %d:%d:%2.3f to %d-%d-%d %d:%d:%2.3f', ...
		s(1),s(2),s(3),s(4),s(5),s(6),e(1),e(2),e(3),e(4),e(5),e(6)))
	end
end
end

if DEBUG, disp(sprintf('\nFound total %d intervals.\n',n_good)), end

if length(var_list)>1
	for j=1:length(var_list)
		eval(['fdm_save{j}=fdm_save_' var_list{j} ';']);
	end
else
	eval(['fdm_save=fdm_save_' var_list{j} ';']);
end
