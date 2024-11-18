function r = log(logLevel,logMsg)
%IRF.LOG   Configurable log routine
%
% IRF.LOG(logLevel) - set the active log level. Default is 'warning'.
%	Enough to specify only the first letter of the log level.
%	logLevel can be 'off','critical','warning','notice' or 'debug'
%
% IRF.LOG(logLevel,logMsg) - output log message logMsg in case the active
%		log level is larger or equal to logLevel. logMsg should be string.
%
% IRF.LOG('log_out',file)     - log output to file.
% IRF.LOG('log_out','screen') - log output to screen, default.
% IRF.LOG('log_out')          - return previously specified file or 'screen'.
% IRF.LOG('short') - output short messages, without debug information
% except when in the debug mode.
% IRF.LOG('full') - output full messages with debug information [default]
%
% Example:
%   irf.log('log_out','/tmp/my_event.log'); % Log to '/tmp/my_event.log'.
%   irf.log('warning'); % set active log level to 'warning'
%                       % prints 'critical' and 'warning' messages
%   irf.log('critical','We should not end in this place of code.')
%   irf.log('warning','Two signals are interpolated')
%   logFile = irf.log('log_out'); % Get name of log file (or 'screen').
%
% See also: LOG_DEBUG_INIT

% Log levels are represented internally by below numbers:
%
% off = 0
% critical = 1
% warning = 2
% notice = 3
% debug = 4

persistent logOut loggingLevel outputShort

if isempty(loggingLevel)
  loggingLevel=2;
end
if isempty(logOut)
  logOut = 'screen';
end
if isempty(outputShort)
  outputShort = false;
end

if nargin == 0
  if nargout
    r = loggingLevel;
  else
    irf.log('warning',log_level_to_msg(loggingLevel,'long'));
  end
  return;
elseif nargin == 1
  if ischar(logLevel)
    switch lower(logLevel(1)) % switch/case is case sensitive
      case 'o'
        loggingLevel = 0;
      case 'c'
        loggingLevel = 1;
      case 'w'
        loggingLevel = 2;
      case 'n'
        loggingLevel = 3;
      case 'd'
        loggingLevel = 4;
      case 's'
        outputShort = true;
      case 'f'
        outputShort = false;
      case 'l'
        if(strcmp(logLevel, 'log_out'))
          % Return log_out (previously specified 'file' or 'screen')
          r = logOut;
          return;
        end
      otherwise
        irf.log('critical','Error! Unrecognized input, see help.');
        error('Unrecognized input.');
    end
    irf.log('warning',['Active log level set to ''' ...
      log_level_to_msg(loggingLevel,'short') '''. '])
    irf.log('warning',log_level_to_msg(loggingLevel,'long'));
  elseif isnumeric(logLevel)
    loggingLevel = logLevel;
    irf.log('warning',log_level_to_msg(loggingLevel,'long'));
  else
    irf.log('critical','Error! Single input parameter should be logLevel, see syntax.');
    error('Wrong syntax');
  end
  return;
end

if loggingLevel==0 % return if level is zero
  return;
end

if nargin == 2
  if ischar(logLevel)
    switch lower(logLevel(1))
      case 'c'
        logLevel = 1;
      case 'w'
        logLevel = 2;
      case 'n'
        logLevel = 3;
      case 'd'
        logLevel = 4;
      case 'l'
        if strcmpi(logLevel,'log_out') % irf.log('log_out',file)
          logOut = logMsg;
          irf.log('warning',['Writing log to ' logOut]);
          return
        end
      otherwise
        irf.log('critical','Error! Unrecognized input, see help.');
        error('Unrecognized input.');
    end
    if logLevel > loggingLevel
      return;
    end
  else
    irf.log('critical','Error! Unrecognized input, see help.');
    error('Unrecognized input.');
  end
else
  irf.log('critical','Error! Max 2 input parameters, see syntax.');
  error('Unrecognized input, max 2 input parameters.');
end

if ~outputShort || ~strcmp(logOut,'screen') || (outputShort && logLevel==4)
  [sta,curr] = dbstack;
  % if irf.log is called from the main env, then use curr,
  % otherwise we are interested in callers name (curr+1)
  if curr == length(sta), idx = curr;
  else, idx = curr +1;
  end
  logMarker = sprintf('%s(%d)',...
    sta(idx).name,...
    sta(idx).line);
  clear sta curr
end

if ~strcmp(logOut,'screen')
  fid = fopen(logOut,'a');
  if fid > 0
    dispStr = ['[' irf_time '][' logMarker '] ' logMsg];
    fprintf(fid,'%s\n',dispStr);
    fclose(fid);
  else
    logOut = 'screen';
    irf.log('critical',['Error! Cannot open output file ' logOut 'for writing'])
    irf.log('critical','Redirecting future output to screen.')
  end
else
  if outputShort && (logLevel ~= 4)
    debugStr = '';
  else
    debugStr = ['[' log_level_to_msg(logLevel,'short') ': ' logMarker '] '];
  end
  disp([repmat(' ',1,(logLevel-1)*2) ...% indentation space depending on level
    debugStr logMsg]);
end
end

function logMsg = log_level_to_msg(logLevel,flag)
% return string describing logLevel
% flag = 'short' (default) 1 > critical
% flag = 'long'            1 > Showing 'critical' messages.
if nargin == 1
  flag = 'short';
end
logMsg = ''; % default
if ischar(flag)
  if strcmpi(flag,'short')
    if logLevel == 1
      logMsg = 'critical';
    elseif logLevel == 2
      logMsg = 'warning';
    elseif logLevel == 3
      logMsg = 'notice';
    elseif logLevel == 4
      logMsg = 'debug';
    end
  elseif strcmpi(flag,'long')
    if logLevel == 1
      logMsg = 'Showing ''critical'' log messages.';
    elseif logLevel == 2
      logMsg = 'Showing ''critical'' and ''warning'' messages.';
    elseif logLevel == 3
      logMsg = 'Showing ''critical'', ''warning'' and ''notice'' messages.';
    elseif logLevel == 4
      logMsg = 'Showing ''critical'', ''warning'', ''notice'' and ''debug'' messages.';
    end
  end
end
end
