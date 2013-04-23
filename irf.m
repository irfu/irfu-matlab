function out=irf(varargin)
% IRF general info on irfu-matlab
%
% IRF runs basic tests
%
% irf('help')
% help irfu-matlab
%		show general help on irfu-matlab.
%
% [out] = irf('check') check if using latest version of irfu-matlab
%	out is logical true if using latest and false if not. 
%
% [out] = irf('mice') check if spice/mice routines are installed properly
% and if necessary add to the path.
% more SPICE info: http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/


if nargin == 0,
	irf('check');
	irf('mice');
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
	case 'demo'
		echodemo irfdemo
	case 'help'
		help irfu-matlab
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
	case 'mice'
		if exist('cspice_j2000','file') % mice is installed
			if (cspice_j2000 == 2451545),
				disp('SPICE/MICE is installed and working properly');
				if nargout, out=true; end
				return;
			else
				disp('SPICE/MICE is installed but NOT WORKING PROPERLY!');
				if nargout, out=false; end
				return;
			end
		else
			disp('adding MICE path to matlab');
			addpath([irf('path') filesep  'mice']);
			ok=irf('mice');
			if ~ok, 
				disp('There are mice problems. Please, contact irfu!');
			end
		end
	case('path')
			out = fileparts(which('irf.m'));	
	case 'check'
		fprintf('Checking if you have latest irfu-matlab...');
		try
			logText      = urlread(logFileUrl);
		catch
			disp('Not connected to internet');
			out = false;
			return;
		end
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
