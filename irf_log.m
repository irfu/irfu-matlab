function irf_log(log_ids,log_msg,varargin)
%IRF_LOG   Configurable logging routine
%
% irf_log(log_ids,log_msg) - log a message LOG_MSG with id LOG_IDS.
% The procedure will check if messages with LOG_IDS are intended for logging
% by current LOG_LEV (see below), and will be send to LOG_OUT (see below).
%
% Input:
% log_id - one of:
%   'fcal' - function calls. Additionally adds function name to the message.
%   'dsrc' - data source (reading of row data from ISDAT and CDF files)
%   'load' - loading of data from MAT files
%   'save' - saving of data to MAT files
%   'calb' - calibration
%   'proc' - processing
% log_msg - message string 
%
% Example:
% irf_log('proc','using E.B=0 for computing Ez')
%
% irf_log('log_lev',level) - set logging level. Default is maximum possible
% log level (all messages are logged).
% LEVEL is a sum of desired options:
%   'fcal' - 1
%   'dsrc' - 2
%   'load' - 4
%   'save' - 8
%   'calb' - 16
%   'proc' - 32
%
% Example:
% irf_log('log_lev',32+16+4)
%   //this allows logging of messages with LOG_IDS 'load','calb' and 'proc'
%
% irf_log('log_out',file) - redirect the output from screen (default) to FILE.
% Default can be restored by running irf_log('log_out','screen')
% 
% Example:
% irf_log('log_out','/tmp/my_event.log')
%
% $Id$

% Copyright 2004,2005 Yuri Khotyaintsev

error(nargchk(2,15,nargin))

%if nargin>2, have_options = 1; args = varargin;
%else, have_options = 0;
%end
%
%while have_options
%	l = 2;
%	if length(args)>=1
%		switch(args{1})
%		case 'sp'
%		
%		otherwise
%        	disp(['Option ''' args{i} '''not recognized'])
%    	end
%		if length(args) > l, args = args(l+1:end);
%		else break
%		end
%	else
%		error('caa:wrongArgType','use c_get_batch(..,''option'',''value'')')
%	end
%end

switch(log_ids)
case 'log_lev'
	global IRF_LOG
	IRF_LOG = log_msg;
	return
case 'log_out'
	global IRF_LOG_OUT
	IRF_LOG_OUT = log_msg;
	return
case 'fcal'
	log_id = 1;
case 'dsrc'
	log_id = 2;
case 'load'
	log_id = 3;
case 'save'
	log_id = 4;
case 'calb'
	log_id = 5;	
case 'proc'
	log_id = 6;
otherwise
	disp('unknown LOG_ID')
	log_id = 255;
end

global IRF_LOG
if isempty(IRF_LOG), log_ok = 1; log_fn = 1;
else
	log_ok = bitget(uint32(IRF_LOG),log_id);
	log_fn = bitget(uint32(IRF_LOG),1);
end

if log_ok
	if log_fn
		[sta,curr] = dbstack;
		log_ids = sprintf('%s(%d) : %s',sta(curr+1).name,sta(curr+1).line,log_ids);
		clear sta curr
	end
	d_str = ['[' datestr(now,31) '][' log_ids '] ' log_msg];
	
	persistent d_out
	persistent d_out_prev
	global IRF_LOG_OUT
	
	if isempty(d_out) | ~strcmp(d_out_prev,IRF_LOG_OUT)
		d_out_prev = IRF_LOG_OUT;
		d_out = 'screen';
		if ~isempty(IRF_LOG_OUT) 
			if ~strcmp(IRF_LOG_OUT,'screen')
				fid = fopen(IRF_LOG_OUT,'a');
				if fid > 0
					d_out = IRF_LOG_OUT;
					fclose(fid);
				else
					disp(['IRF_LOG: cannot open output file ' IRF_LOG_OUT 'for writing'])
					disp('IRF_LOG: check the global variable IRF_LOG_OUT')
					disp('IRF_LOG: redirecting output to screen')
				end
			end
		end
	end
	if strcmp(d_out,'screen')
		disp(d_str)
	else
		fid = fopen(d_out,'a');
		if fid > 0
			fprintf(fid,'%s\n',d_str);
			fclose(fid);
		else 
			disp(['IRF_LOG: cannot open output file ' IRF_LOG_OUT 'for writing'])
			disp('IRF_LOG: redirecting output to screen')
			disp(d_str)
		end
	end
end
