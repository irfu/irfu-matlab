function c_log(log_ids,log_msg,varargin)
%C_LOG configurable logging routine
%
% c_log(log_ids,log_msg) - log a message LOG_MSG with id LOG_IDS.
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
% c_log('proc','using E.B=0 for computing Ez')
%
% c_log('log_lev',level) - set logging level. Default is maximum possible
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
% c_log('log_lev',32+16+4)
%   //this allows logging of messages with LOG_IDS 'load','calb' and 'proc'
%
% c_log('log_out',file) - redirect the output from screen (default) to FILE.
% Default can be restored by running c_log('log_out','screen')
% 
% Example:
% c_log('log_out','/tmp/my_event.log')
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev

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
%		error('caa:wrongArgType','use getRData(..,''option'',''value'')')
%	end
%end

switch(log_ids)
case 'log_lev'
	global C_LOG
	C_LOG = log_msg;
	return
case 'log_out'
	global C_LOG_OUT
	C_LOG_OUT = log_msg;
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

global C_LOG
if isempty(C_LOG), log_ok = 1; log_fn = 1;
else
	log_ok = bitget(uint32(C_LOG),log_id);
	log_fn = bitget(uint32(C_LOG),1);
end

if log_ok
	if log_fn, log_ids = [evalin('caller','mfilename') ' : ' log_ids]; end
	d_str = ['[' datestr(now,31) '][' log_ids '] ' log_msg];
	
	persistent d_out
	persistent d_out_prev
	global C_LOG_OUT
	
	if isempty(d_out) | ~strcmp(d_out_prev,C_LOG_OUT)
		d_out_prev = C_LOG_OUT;
		d_out = 'screen';
		if ~isempty(C_LOG_OUT) 
			if ~strcmp(C_LOG_OUT,'screen')
				fid = fopen(C_LOG_OUT,'a');
				if fid > 0
					d_out = C_LOG_OUT;
					fclose(fid);
				else
					disp(['C_LOG: cannot open output file ' C_LOG_OUT 'for writing'])
					disp('C_LOG: check the global variable C_LOG_OUT')
					disp('C_LOG: redirecting output to screen')
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
			disp(['C_LOG: cannot open output file ' C_LOG_OUT 'for writing'])
			disp('C_LOG: redirecting output to screen')
			disp(d_str)
		end
	end
end
