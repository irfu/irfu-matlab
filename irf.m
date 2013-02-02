function out=irf(varargin)
% IRF general info on irfu-matlab
%
% to obtain help on IRFU matlab routines run "irf" or
%  >> help irfu-matlab
if nargin == 0,
	irf('check');
	help irfu-matlab;
	return;
else
	if ischar(varargin{1}),
		action = lower(varargin{1});
	else
		irf_log('fcal','unknown input parameters');
	end
end

logFileUrl = 'http://www.space.irfu.se/~andris/irfu-matlab/log.txt';

switch action
	case 'version'
		fid=fopen(log_file);
        tline = fgetl(fid);
		fclose(fid);
		versionTime = tline(1:10);
		if nargout == 0,
			disp(['You are using irfu-matlab version from ' versionTime]);
		else
			out = versionTime;
		end
	case 'check'
		fprintf('Checking if you have latest irfu-matlab...');
		logText      = urlread(logFileUrl);
		logTextArray = textscan(logText, '%s', 'delimiter', sprintf('\n'));
		logTextArray = logTextArray{1};
		newestVersion = logTextArray{1}(1:10);
		currentVersion = irf('version');
		if ~strcmp(newestVersion,currentVersion)
			indices = find(cellfun(@(x) any(strcmp(x(1:10),currentVersion)),logTextArray));
			if indices > 1,
				disp('NO!');
				disp(' ');
				disp(['Newest irfu-matlab is from: ' newestVersion]);
				disp(['  Your irfu-matlab is from: ' currentVersion]);
				disp('Please update, see <a href="https://launchpad.net/irfu-matlab">https://launchpad.net/irfu-matlab</a>');
				disp('Log of updates: ');
				disp(logTextArray(1:indices-1));
				disp(' ');
			end
			if nargout, out = false; end 
		else
			disp('YES:)');
			if nargout, out = true; end
		end
	otherwise
		irf_log('fcal','unknown argument');
end

function logFile=log_file
irfDir = fileparts(which('irf_plot.m'));
logFile = [irfDir filesep 'log.txt'];
