function out=irf(varargin)
% IRF general info on irfu-matlab
%
% IRF checks version, installed libraries and sets up necessary pathes
%
% IRF('help') or 'help irfu-matlab' shows general help on irfu-matlab.
%
% [out] = IRF('check') check if using latest version of irfu-matlab
%	out is logical true if using latest and false if not.
%
% [out] = IRF('mice') check if spice/mice routines are installed properly
% and if necessary add to the path. run irf('mice_help') if you want to see
% more help on mice kernels.
% more SPICE info: http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/
%
% [out] = IRF('irbem')
%
% Check if ONERA IRBEM library is installed
% http://craterre.onecert.fr/prbem/irbem/description.html
%
% [out] = IRF('ceflib')
%
% Check if IRAP CEFLIB is installed, http://ceflib.irap.omp.eu/
%
% version = IRF('version') return IRF version
%
% IRF('demo') demonstration how to use IRF

%% Defaults
% file to check version
logFileUrl = 'https://raw.github.com/irfu/irfu-matlab/master/log.txt';

%% Input check
if nargin == 0,
	irf('check');
	irf('ceflib');
	irf('mice');
	irf('irbem');
	irf('check_path');
	return;
else
	if ischar(varargin{1}),
		action = lower(varargin{1});
	else
		irf_log('fcal','unknown input parameters');
	end
end

%% Actions
switch lower(action)
	case 'check'
		fprintf('Checking if you have latest irfu-matlab... ');
		try
			logText      = urlread(logFileUrl);
		catch
			disp('Not connected to internet');
			out = false;
			return;
		end
		logTextArray = textscan(logText, '%s', 'delimiter', sprintf('\n'));
		logTextArray = logTextArray{1};
		iSpace = strfind(logTextArray{1},' ');
		newestVersion = logTextArray{1}(iSpace(1):iSpace(2)-1);
		currentVersion = irf('version');
		if ~strcmp(newestVersion,currentVersion)
			indices = find(cellfun(@(x) any(strfind(x,currentVersion)),logTextArray));
			if indices > 1,
				disp('NO!');
				disp(' ');
				disp(['Newest irfu-matlab is from: ' newestVersion]);
				disp(['  Your irfu-matlab is from: ' currentVersion]);
				disp('Please update, see <a href="https://github.com/irfu/irfu-matlab">https://github.com/irfu/irfu-matlab</a>');
				disp('Log of updates: ');
				for iInd = 1 : indices -1
					fprintf('%s\n',logTextArray{iInd})
				end
				disp(' ');
      else
        disp('You are at the bleeding edge :-)');
			end
			if nargout, out = false; end
		else
			disp('YES:)');
			if nargout, out = true; end
		end
	case 'check_path'
		irfPath = [irf('path') filesep];
		notOnIrfPath = @(x) ~any(strfind(path, [irfPath x]));
		contribDirectories = {...
			['contrib' filesep 'isdat'],...
			['contrib' filesep 'libirbem'],...
			['contrib' filesep 'libcef'],...
			['contrib' filesep 'matlab_central'],...
			['contrib' filesep 'matlab_central' filesep 'cm_and_cb_utilities'],...
			['contrib' filesep 'mice'],...
      ['contrib' filesep 'nasa_cdf_patch'],...
			};
		irfDirectories = {'irf','plots',...
			['mission' filesep 'cluster'],...
			['mission' filesep 'solar_orbiter'],...
			['mission' filesep 'cluster' filesep 'caa'],...
			};
		strPath = {contribDirectories{:},irfDirectories{:}}; %#ok<CCAT>
		for iPath = 1:numel(strPath)
			if notOnIrfPath(strPath{iPath}),
				pathToAdd = [irfPath strPath{iPath}];
				addpath(pathToAdd);
				disp(['Added to path: ' pathToAdd]);
			end
		end
	case 'demo'
		echodemo irfdemo
	case 'help'
		help irfu-matlab
	case 'mice'
		if exist('cspice_j2000','file') % mice is installed
			if (cspice_j2000 == 2451545),
				disp('SPICE/MICE is OK');
				if nargout, out=true; end
				return;
			else
				disp('SPICE/MICE is installed but NOT WORKING PROPERLY!');
				if nargout, out=false; end
				return;
			end
		else
			micePath = [irf('path') filesep 'contrib' filesep  'mice'];
			disp(['adding MICE path to matlab: ' micePath]);
			addpath(micePath);
			ok=irf('mice');
			if ~ok,
				disp('There are mice problems. Please, contact irfu!');
			end
		end
	case 'mice_help'
		disp('Kernel files at IRFU are located at spis:/share/SPICE');
		disp('Kernels at irfu: general, Cassini, Rosetta, Solar Orbiter, JUICE');
		disp('');
		disp('At other locations if you want to get kernel files, create directory and run:');
		disp('> wget  --timestamping -r -nH --cut-dirs=2 -X *a_old_versions* ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels');
		disp('This will download all the latest versions of necessary general kernels.');
		disp('If you want for example get all Rosetta kernels, execute:');
		disp('> wget  --timestamping -r -nH --cut-dirs=2 -X *former_versions* ftp://naif.jpl.nasa.gov/pub/naif/ROSETTA');
		disp('');
	case 'irbem'
		if exist('onera_desp_lib_coord_trans','file') % irbem is installed
			x=[0 0 1];
			try
				y=onera_desp_lib_coord_trans([0 0 1],'gse2geo', now);
				yy=onera_desp_lib_coord_trans(y,'geo2gse',now);
			catch
				disp('IRBEM is not properly installed!');
				if nargout, out=false; end
				return;
			end
			if (max(abs(yy-x))<1e-3),
				disp('IRBEM is OK');
				if nargout, out=true; end
				return;
			else
				disp('IRBEM is installed but NOT WORKING PROPERLY!');
				disp('gse>geo>gse differs by more than 0.1% from original vector');
				if nargout, out=false; end
				return;
			end
		else
			oneraPath = [irf('path') filesep 'contrib' filesep  'libirbem'];
			disp(['adding IRBEM path to matlab: ' oneraPath]);
			addpath(oneraPath);
			ok=irf('irbem');
			if ~ok,
				disp('There are IRBEM problems. Please, contact irfu!');
			end
		end
	case 'ceflib'
		if exist('cef_init','file') % CESR CEFLIB is installed
			cef_init();
			cef_verbosity(0);
			if ( cef_read(which('C1_CP_EFW_L3_P__20010201_120000_20010201_120100_V110503.cef.gz'))==0 && ...
					numel(cef_date(cef_var ('time_tags'))) == 15 && ...
					numel(cef_var('Spacecraft_potential')) == 15 )
				disp('CEFLIB is OK');
				if nargout, out = true; end
			else
				disp('There are CEFLIB problems. Please, contact irfu!');
				if nargout, out = false; end
			end
		else
			ceflibPath = [irf('path') filesep 'contrib' filesep  'libcef'];
			disp(['adding CEFLIB path to matlab: ' ceflibPath]);
			addpath(ceflibPath);
			out=irf('ceflib');
			if ~out,
				disp('There are CEFLIB problems. Please, contact irfu!');
			end
		end
	case 'path'
		out = fileparts(which('irf.m'));
	case 'version'
		logFile = [fileparts(which('irf.m')) filesep 'log.txt'];
		fid=fopen(logFile);
		tline = fgetl(fid);
		fclose(fid);
		iSpace = strfind(tline,' ');
		versionTime   = tline(1:iSpace(1)-1);
		versionNumber = tline(iSpace(1):iSpace(2)-1);
		if nargout == 0,
			disp(['irfu-matlab version: ' versionTime ', ' versionNumber]);
		else
			out = versionNumber; % return only date
		end
	otherwise
		irf_log('fcal','unknown argument');
end


