function r = log(logLevel,logMsg)
%IRF.LOG   Configurable logging routine
%
% IRF.LOG(level) - set logging level. Default is 1, if zero do not show log.
% 
% IRF.LOG(logLevel,logMsg) - log message logMsg if current
%		logging level is larger or equal to logLevel.
%
% logLevel - can have values 1,2,..
% logMsg   - message string
%
% IRF.LOG(logMsg) - log message logMsg with logLevel=1
% 
% Example:
%   irf.log(2); % print level 1 & 2 messages
%   irf.log(1,'Very important message')
%   irf.log(2,'Less important message')
%
% IRF.LOG('log_out',file)     - log output to file.
% IRF.LOG('log_out','screen') - log output to screen, default. 
%
% Example:
%	irf.log('log_out','/tmp/my_event.log')

persistent logOut
persistent loggingLevel
if isempty(loggingLevel),
    loggingLevel=1;
end
if isempty(logOut)
	logOut = 'screen';
end

if nargin == 0,
  if nargout, r = loggingLevel;
  else
    irf.log(1,['Current logging level is ' num2str(loggingLevel)]);
  end
	return;
elseif nargin == 1, 
	if isnumeric(logLevel),
		loggingLevel = logLevel;
		irf.log(2,['logging level set to ' num2str(logLevel)]);
	elseif ischar(logLevel) % irf.log(1,logMsg)
		irf.log(1,logLevel);
		return;
	else
		irf.log(1,'Warning! Single input parameter should be number, see syntax.');
	end
	return;
end

if loggingLevel==0 % return if level is zero
	return;
end

if nargin == 2
	if isnumeric(logLevel),  % irf.log(logLevel,logMsg)
		if logLevel > loggingLevel,
			return;
		end
	elseif ischar(logLevel) && strcmpi(logLevel,'log_out'), % irf.log('log_out',file)
		logOut = logMsg;
		irf.log(2,['Writing log to ' logOut]);
		return
	else
		irf.log(1,'Warning! Unrecognized input, see help.');
		return;
	end
else
	irf.log(1,'Warning! Max 2 input parameters, see syntax.');
	return;
end	

[sta,curr] = dbstack;
% if irf.log is called from the main env, then use curr,
% otherwise we are interested in callers name (curr+1)
if curr == length(sta), idx = curr;
else idx = curr +1;
end
logMarker = sprintf('%s(%d) : %d',...
	sta(idx).name,...
	sta(idx).line,...
	logLevel);
clear sta curr

if ~strcmp(logOut,'screen')
	fid = fopen(logOut,'a');
	if fid > 0
		dispStr = ['[' irf_time '][' logMarker '] ' logMsg];
		fprintf(fid,'%s\n',dispStr);
		fclose(fid);
	else
		logOut = 'screen';
		irf.log(1,['Error! Cannot open output file ' logOut 'for writing'])
		irf.log(1,'Redirecting future output to screen.')
	end
else
	dispStr = [repmat(' ',1,(logLevel-1)*2) ...% indentation space depending on level
		'[' num2str(logLevel) '] ' logMsg];
	disp(dispStr)
end
