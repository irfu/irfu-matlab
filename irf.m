function [out,out1]=irf(varargin)
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
% more SPICE info: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/
%
% [out] = IRF('irbem')
%
% Check if ONERA IRBEM library is installed
% https://craterre.onecert.fr/prbem/irbem/description.html
%
% [out] = IRF('ceflib')
%
% Check if IRAP CEFLIB is installed, http://ceflib.irap.omp.eu/
%
% version = IRF('version') return IRF version number
% [versionNumber, versionDate] = IRF('version') return also date
%
% [out] = IRF('path') returns the path to the current irfu-matlab root
% directory.
%
% IRF('demo') demonstration how to use IRF

%this is an edit to load BLAS
ones(10)*ones(10); %#ok<VUNUS>

%% Defaults
% file to check version
logFileUrl = 'https://raw.githubusercontent.com/irfu/irfu-matlab/master/log.txt';

%% Input check
if nargin == 0
  setenv('LC_ALL','C'); % temporar fix for R2014a problems on Unix, https://www.mathworks.com/matlabcentral/answers/126994-locale-system-error
  irf('check_path');
  irf('check');
  irf('ceflib');
  irf('mice');
  irf('irbem');
  irf('check_os');
  irf('matlab');
  irf('cdf_leapsecondstable');
  return;
else
  if ischar(varargin{1})
    action = lower(varargin{1});
  else
    irf.log('critical','string input required');
    error('string input required')
  end
end

%% Actions
switch lower(action)
  case 'check'
    [currentVersion,currentVersionDate] = irf('version');
    disp(['irfu-matlab version: ' currentVersion]);
    fprintf('Checking if you have latest irfu-matlab... ');
    try
      if verLessThan('matlab', '8.4')
        logText = urlread(logFileUrl); %#ok<URLRD> webread introduced in R2014b
      else
        logText = webread(logFileUrl);
      end
    catch ME
      disp(['Failed to get upstream version information, resulted in error message: ', ME.message]);
      disp(['  Your irfu-matlab: ' currentVersion ...
        ' from ' currentVersionDate]);
      out = false;
      return;
    end
    logTextArray = textscan(logText, '%s', 'delimiter', sprintf('\n')); %#ok<SPRINTFN>
    logTextArray = logTextArray{1};
    iSpace = strfind(logTextArray{1},' ');
    newestVersion = logTextArray{1}(iSpace(1):iSpace(2)-1);
    if ~strcmp(newestVersion,currentVersion)
      indices = find(cellfun(@(x) any(strfind(x,currentVersion)),logTextArray));
      if indices > 1
        disp('NO!');
        disp(' ');
        disp(['Newest irfu-matlab version: ' newestVersion]);
        disp('Please update, see <a href="https://github.com/irfu/irfu-matlab">https://github.com/irfu/irfu-matlab</a>');
        disp('Log of updates: ');
        for iInd = 1 : indices -1
          fprintf('%s\n',logTextArray{iInd})
        end
        disp(' ');
      else
        disp('unclear! you are either the bleeding edge or have to update :-)');
        disp(['Newest irfu-matlab is from: ' newestVersion]);
        disp(['  Your irfu-matlab is from: ' currentVersion]);
      end
      if nargout, out = false; end
    else
      disp('YES :)');
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
      ['plots'   filesep 'mms'],...
      ['plots'   filesep 'solo'],...
      ['mission' filesep 'cluster'],...
      ['mission' filesep 'cluster' filesep 'caa'],...
      ['mission' filesep 'juice'],...
      ['mission' filesep 'solar_orbiter'],...
      ['mission' filesep 'solar_orbiter' filesep 'bicas' filesep 'src'],...
      ['mission' filesep 'themis'],...
      ['mission' filesep 'thor'],...
      ['mission' filesep 'mms'],...
      ['mission' filesep 'mms' filesep 'mms_testFunctions'],...
      };
    strPath = {contribDirectories{:},irfDirectories{:}}; %#ok<CCAT>
    if ~any(strfind(path, irf('path'))) % irfu-matlab root folder
      addpath(irfPath);
      disp(['Added to path: ' irfPath]);
    end
    for iPath = 1:numel(strPath)
      if notOnIrfPath(strPath{iPath})
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
      try
        if(cspice_j2000 == 2451545)
          disp('SPICE/MICE is OK');
          if nargout, out = true; end
          return
        else
          disp('SPICE/MICE is installed but NOT WORKING PROPERLY!');
          if nargout, out = false; end
          return
        end
      catch
        disp('SPICE/MICE is installed but NOT WORKING PROPERLY!');
        if nargout, out = false; end
        return
      end
    else
      micePath = [irf('path') filesep 'contrib' filesep  'mice'];
      disp(['adding MICE path to matlab: ' micePath]);
      addpath(micePath);
      ok = irf('mice');
      if ~ok
        disp('MICE  .. NOT OK. Please contact IRFU if you need SPICE/MICE for your intended use of irfu-matlab!');
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
        t0 = now;
        y=onera_desp_lib_coord_trans([0 0 1],'gse2geo', t0);
        yy=onera_desp_lib_coord_trans(y,'geo2gse',t0);
        if (max(abs(yy-x))<1e-3)
          disp('IRBEM is OK');
          if nargout, out=true; end
          return;
        else
          disp('IRBEM is installed but NOT WORKING PROPERLY!');
          disp('gse>geo>gse differs by more than 0.1% from original vector');
          if nargout, out=false; end
          return;
        end
      catch
        if ismac
          disp('IRBEM .. not OK. Please check that:');
          disp('1) you have Xcode installed');
          disp('2) you have GFortran installed');
          disp('   if not, download GF 5.1 from http://hpc.sourceforge.net');
          disp('3) open /Applications/MATLAB_R2017b.app/bin/mexopts.sh');
          disp('   check that in maci64 section MAC OS version number');
          disp('   corresponds to your OS version, e.g. 10.9 for Mavericks.');
          disp('   run in matlab > mex -setup');
          disp('Or see an alternative workaround here:');
          disp('  <a href="https://github.com/irfu/irfu-matlab/issues/61#issuecomment-797499907">https://github.com/irfu/irfu-matlab/issues/61#issuecomment-797499907</a>');
        elseif ispc
          disp('IRBEM .. not OK. If this package is required for your intended use of irfu-matlab:');
          disp('   Please follow the installation instructions for Windows');
          disp('   https://sourceforge.net/p/irbem/code/HEAD/tree/trunk/manual/frames/matlab.html?format=raw');
          disp('   and install the package and all required libraries into: ');
          disp(['   ',irf('path'),filesep,'contrib',filesep,'libirbem']);
        else
          disp('IRBEM .. not OK. Please contact IRFU if you need IRBEM for your intended use of irfu-matlab!');
        end
        if nargout, out=false; end
        return;
      end
    else
      oneraPath = [irf('path') filesep 'contrib' filesep  'libirbem'];
      disp(['adding IRBEM path to matlab: ' oneraPath]);
      addpath(oneraPath);
      ok=irf('irbem');
      if ~ok
        disp('IRBEM .. NOT OK. Please contact IRFU if you need IRBEM for your intended use of irfu-matlab!');
      end
    end

  case 'ceflib'
    if ~ispc
      if exist('cef_init','file') % CESR CEFLIB is installed
        try % Try but don't crash on failure (see irfu-matlab issue #66)
          cef_init();
          cef_verbosity(0);
          if ( cef_read(which('C1_CP_EFW_L3_P__20010201_120000_20010201_120100_V110503.cef.gz'))==0 && ...
              numel(cef_date(cef_var ('time_tags'))) == 15 && ...
              numel(cef_var('Spacecraft_potential')) == 15 )
            disp('CEFLIB is OK');
            if nargout, out = true; end
            datastore('irfu_matlab','okCeflib',true);
          else
            disp('There are CEFLIB problems. Please contact IRFU if you need CEF for your intended use of irfu-matlab!');
            if nargout, out = false; end
            datastore('irfu_matlab','okCeflib',false);
          end
        catch
          disp('There are CEFLIB problems. Please contact IRFU if you need CEF for your intended use of irfu-matlab!');
          if nargout, out = false; end
          datastore('irfu_matlab','okCeflib',false);
        end
      else
        ceflibPath = [irf('path') filesep 'contrib' filesep  'libcef'];
        disp(['adding CEFLIB path to matlab: ' ceflibPath]);
        addpath(ceflibPath);
        out=irf('ceflib');
        if ~out
          disp('There are CEFLIB problems. Please, contact irfu!');
          datastore('irfu_matlab','okCeflib',false);
        end
      end
    else
      datastore('irfu_matlab','okCeflib',false);
    end

  case 'cdf_leapsecondstable'
    % Check to see if CDF_LEAPSECONDSTABLE is set as environment
    % variable, as it is used by TT2000 conversions in irf_time. If it is
    % not, try to set it to the included CDFLeapSeconds.txt. If it is set
    % verify it points to an exisiting file which is up to date with the
    % included CDFLeapSeconds.txt. If the file it points to does
    % not exist try to set it to the included file instead.
    leapsecondstable = getenv('CDF_LEAPSECONDSTABLE');
    if( isempty(leapsecondstable) )
      % It was not set, try to set it.
      disp('CDF_LEAPSECONDSTABLE was not set in user environment.');
      disp(['Automatically setting it to: ',which('CDFLeapSeconds.txt')]);
      if(nargout), out=true; end
      datastore('irfu_matlab','okLeapsecondstable',true);
      try setenv('CDF_LEAPSECONDSTABLE',which('CDFLeapSeconds.txt'));
      catch
        disp('Was unsuccesful in setting environment variable.');
        disp('Falling back to hard coded leap seconds table.');
        if(nargout), out=false; end
        datastore('irfu_matlab','okLeapsecondstable',false);
      end
    else
      % It was set, check to see if it is a valid file and if it is up to
      % date with version in irfu-matlab/contrib/nasa_cdf_patch.
      if( exist(leapsecondstable,'file') )
        fileId = fopen(which('CDFLeapSeconds.txt'),'r');
        leapIRFU = textscan(fileId,'%s','delimiter','\n');
        fclose(fileId);
        leapIRFU = leapIRFU{1,1}(end); % Line of last leap second added
        leapIRFU = str2double(strsplit(leapIRFU{1},' ')); % Convert to num
        fileId = fopen(leapsecondstable,'r');
        leapLocal = textscan(fileId,'%s','delimiter','\n');
        fclose(fileId);
        leapLocal = leapLocal{1,1}(end); % Line of last leap second added
        leapLocal = str2double(strsplit(leapLocal{1},' ')); % Convert to num
        if(leapLocal(4)<leapIRFU(4))
          disp('Local leap seconds table appears to be out of date.');
          disp('Please check you system settings and update your table.');
          disp(['Your environment variable CDF_LEAPSECONDSTABLE is set to ',leapsecondstable]);
          fprintf('Presently your table contain leap second %d-%d-%d.\n',leapLocal(1),leapLocal(2),leapLocal(3));
          fprintf('while the latest leap second is %d-%d-%d.\n',leapIRFU(1),leapIRFU(2),leapIRFU(3));
          disp('For now, we will use your specified leap seconds table.');
          if(nargout), out=true; end
          datastore('irfu_matlab','okLeapsecondstable',true);
        else
          disp('CDF_LEAPSECONDSTABLE is OK');
          if(nargout), out=true; end
          datastore('irfu_matlab','okLeapsecondstable',true);
        end
      else
        disp(['Your environment variable CDF_LEAPSECONDSTABLE was set to ',leapsecondstable,' but this file does not exist?']);
        disp(['Automatically setting it to: ',which('CDFLeapSeconds.txt'),' instead.']);
        if(nargout), out=true; end
        datastore('irfu_matlab','okLeapsecondstable',true);
        try setenv('CDF_LEAPSECONDSTABLE',which('CDFLeapSeconds.txt'));
        catch
          disp('Was unsuccesful in setting environemnt variable, falling back to hard coded leap seconds table.');
          if(nargout), out=false; end
          datastore('irfu_matlab','okLeapsecondstable',false);
        end
      end
    end

  case 'check_os'
    % NASA's SPDF cdf patch use compiled mex files. For irfu-matlab only
    % Linux, Mac and Windows (all of which 64 bit) OS are included.
    switch(computer)
      case({'GLNXA64','PCWIN64','MACI64'})
        % OK, system is supported.
        disp('Operating system is OK');
        if(nargout), out=true; end
        datastore('irfu_matlab','okCheckOS',true);
      case('PCWIN')
        % Untested, no guarantee it will work.
        disp('Operating system, Windows 32 bit, is not fully tested with IRFU-MATLAB.');
        disp('Recommended OS are: Linux, Mac and Windows (all 64 bit).');
        if(nargout), out=true; end
        datastore('irfu_matlab','okCheckOS',true);
      otherwise
        disp('Currently only compiled SPDF* mex files for the Linux, Windows and Mac operating systems (all 64bit) are included.');
        disp('If running on other system, please contact IRFU for help.');
        if(nargout), out=false; end
        datastore('irfu_matlab','okCheckOS',false);
    end

  case 'matlab'
    % Issue warning if running too old Matlab. This should be incremented
    % when irfu-matlab relies on newer Matlab functions not found in older
    % versions of Matlab.
    if(verLessThan('matlab','8.4'))
      warning(['IRFU-Matlab relies on code introduced in Matlab R2014b, ',...
        'please look into upgrading your Matlab installation or contacting IRFU for help.']);
    else
      disp('Matlab version is OK');
      if(nargout), out=true; end
      datastore('irfu_matlab','okMatlab',true);
    end

  case 'path'
    out = fileparts(which('irf.m'));

  case 'version'
    logFile = [fileparts(which('irf.m')) filesep 'log.txt'];
    fid = fopen(logFile);
    tline = fgetl(fid);
    fclose(fid);
    iSpace = strfind(tline,' ');
    versionTime   = tline(1:iSpace(1)-1);
    versionNumber = tline(iSpace(1):iSpace(2)-1);
    if nargout == 0
      disp(['irfu-matlab version: ' versionTime ', ' versionNumber]);
    else
      out = versionNumber; % return only date
      out1 = versionTime;
    end

  otherwise
    error('unknown input argument');
end

