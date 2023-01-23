function irf_log(log_ids,log_msg)
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
%   IRF_LOG ON  - swith ON the showing of log messages
%   IRF_LOG OFF - switch OFF the showing of log messages
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

persistent encourageSwitchToIrfLog
persistent d_out
persistent d_out_prev
persistent IRF_LOG IRF_LOG_OUT
persistent flag_irf_log_ON
if isempty(encourageSwitchToIrfLog)
  disp('  ');
  disp(' IRF.LOG prefered instead of IRF_LOG');
  disp(' Syntax is slightly changed, see help');
  disp('  ');
  encourageSwitchToIrfLog = false;
end

if isempty(flag_irf_log_ON)
  flag_irf_log_ON=1;
end
if ~flag_irf_log_ON % if irf_log is off return
  return;
end

narginchk(1,15)

switch lower(log_ids)
  case 'log_lev'
    IRF_LOG = log_msg;
    return
  case 'log_out'
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
  case 'on'
    flag_irf_log_ON=1;
    return
  case 'off'
    flag_irf_log_ON=0;
    return
  otherwise
    [sta,curr] = dbstack;
    % if irf_log is called from the main env, then use curr,
    % otherwise we are interested in callers name (curr+1)
    if curr == length(sta), idx = curr;
    else, idx = curr +1;
    end
    log_ids = sprintf('%s(%d) : %s',sta(idx).name,sta(idx).line,log_ids);
    disp(['unknown LOG_ID at ' log_ids])
    log_id = 255;
end

if isempty(IRF_LOG), log_ok = 1; log_fn = 1;
else
  log_ok = bitget(uint32(IRF_LOG),log_id);
  log_fn = bitget(uint32(IRF_LOG),1);
end

if log_ok
  if log_fn
    [sta,curr] = dbstack;
    % if irf_log is called from the main env, then use curr,
    % otherwise we are interested in callers name (curr+1)
    if curr == length(sta), idx = curr;
    else, idx = curr +1;
    end
    log_ids = sprintf('%s(%d) : %s',sta(idx).name,sta(idx).line,log_ids);
    clear sta curr
  end
  
  if isempty(d_out) || ~strcmp(d_out_prev,IRF_LOG_OUT)
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
  if flag_irf_log_ON
    if strcmp(d_out,'screen') % write to screen, do not include time
      d_str = ['[' log_ids '] ' log_msg];
      disp(d_str)
    else % write to file, include time in log string
      d_str = ['[' irf_time '][' log_ids '] ' log_msg];
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
end
