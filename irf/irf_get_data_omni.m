function f = irf_get_data_omni( tint, parameter , database)
%IRF_GET_DATA_OMNI download omni data
%
% f=IRF_GET_DATA_OMNI(tint,parameter) download parameters for specified time interval
%       f is column vector, first column time and parameters in next columns
%
% parameters - string, paramters separated by comma (case does not matter).
%               'B'     - <|B|>, magnetic field magnitude [nT]
%               'avgB'  - |<B>|, magnitude of average magnetic field [nT]
%               'Blat'  - latitude angle of average B field (GSE)
%               'Blong' - longitude angle of average B field (GSE)
%               'Bx'    - Bx GSE (the same as 'BxGSE')
%               'By'
%               'Bz'
%               'ByGSM' - By GSM
%               'BzGSM' - Bz GSM
%               'T'     - proton temperature (K)
%               'n'     - proton density (cc)
%               'NaNp'  - alpha/proton ratio
%               'v'     - bulk speed (km/s)
%               'vx'    - vx GSE
%               'vy'    - vy GSE, corrected for abberation (29.8 km/s)
%               'vz'    - vz GSE
%               'vlon'  - bulk flow longitude in degrees
%               'vlat'  - bulk flow latitude in degrees
%               'E'     - electric field (mV/m)
%               'P'     - flow pressure (nPa)
%               'beta'  - plasma beta
%               'Ma'    - Alfven Mach number
%               'bsnx','bsny','bsnz' - x,y,z position of bow shock nose
%               'Ms'    - 1 AU IP Magnetosonic Mach number
%               'ssn'   - daily sunspot number
%               'dst'   - DST index
%               'kp'    - Kp index
%               'f10.7' - F10.7 flux
%               'pc'    - PC index
%               'ae'    - AE index
%               'al'    - AL index
%               'au'    - AL index
%               'imfid' - Spacecraft ID for IMF, 50=IMP8, 51=WIND,
%                         60=Geotail, 71=ACE
%               'swid'  - Spacecraft ID for SW
%               'ts'    - Timeshift to bow shock in seconds
%               'rmsts' - Root mean square of timeshift to bow shock
%
% f=IRF_GET_DATA_OMNI(tint,parameter,database) download from specified database
%
% database:  'omni_hour' or 'omni2' - 1h resolution OMNI data (default)
%            'omni_min'             - 1min resolution OMNI data
%
% Examples:
%   tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 12 31 23 59 0])];
%   ff= irf_get_data_omni(tint,'b,bx,bygsm');
%   ff= irf_get_data_omni(tint,'f10.7');
%   ff= irf_get_data_omni(tint,'b','omni_min');

% https://omniweb.gsfc.nasa.gov/html/ow_data.html

% data description of high res 1min,5min omni data set
% https://omniweb.gsfc.nasa.gov/html/omni_min_data.html#4b
% %
% %
% Year                        I4 1995 ... 2006
% Day                         I4 1 ... 365 or 366
% Hour                        I3 0 ... 23
% Minute                      I3 0 ... 59 at start of average
% ID for IMF spacecraft	      I3 See  footnote D below
% ID for SW Plasma spacecraft	I3 See  footnote D below
% # of points in IMF averages	I4
% # of points in Plasma averages	I4
% Percent interp		          I4 See  footnote A above
% Timeshift, sec		          I7
% RMS, Timeshift		          I7
% RMS, Phase front normal	    F6.2 See Footnotes E, F below
% Time btwn observations, sec	I7 DBOT1, See  footnote C above
% Field magnitude average, nT	F8.2
% Bx, nT (GSE, GSM)		        F8.2
% By, nT (GSE)		            F8.2
% Bz, nT (GSE)		            F8.2
% By, nT (GSM)	              F8.2 Determined from post-shift GSE components
% Bz, nT (GSM)	              F8.2 Determined from post-shift GSE components
% RMS SD B scalar, nT	        F8.2
% RMS SD field vector, nT	    F8.2 See  footnote E below
% Flow speed, km/s		        F8.1
% Vx Velocity, km/s, GSE	    F8.1
% Vy Velocity, km/s, GSE      F8.1
% Vz Velocity, km/s, GSE      F8.1
% Proton Density, n/cc    		F7.2
% Temperature, K		          F9.0
% Flow pressure, nPa	       	F6.2 See  footnote G below
% Electric field, mV/m		    F7.2 See  footnote G below
% Plasma beta		              F7.2 See  footnote G below
% Alfven mach number          F6.1 See  footnote G below
% X(s/c), GSE, Re		          F8.2
% Y(s/c), GSE, Re		          F8.2
% Z(s/c), GSE, Re		          F8.2
% BSN location, Xgse, Re      F8.2 BSN = bow shock nose
% BSN location, Ygse, Re	    F8.2
% BSN location, Zgse, Re 	    F8.2
%
% AE-index, nT                I6      See World Data Center for Geomagnetism, Kyoto
% AL-index, nT                I6      See World Data Center for Geomagnetism, Kyoto
% AU-index, nT                I6      See World Data Center for Geomagnetism, Kyoto
% SYM/D index, nT             I6      See World Data Center for Geomagnetism, Kyoto
% SYM/H index, nT             I6      See World Data Center for Geomagnetism, Kyoto
% ASY/D index, nT             I6      See World Data Center for Geomagnetism, Kyoto
% ASY/H index, nT             I6      See World Data Center for Geomagnetism, Kyoto
% PC(N) index,                F7.2    See World Data Center for Geomagnetism, Copenhagen
% Magnetosonic mach number    F5.1    See  footnote G below

%% Define selected database
if nargin < 3 % database not specified defaulting to omni2
  omniDatabase    = 'omni_hour';
elseif nargin == 3 && ischar(database) % database specified
  if strcmpi(database,'omni2') || strcmpi(database,'omni_hour')
    omniDatabase  = 'omni_hour';
  elseif strcmpi(database,'omni_min') || strcmpi(database,'min')
    omniDatabase  ='omni_min';
  else
    errStr = ['Unknown database: ' database];
    irf.log('critical',errStr); error(errStr);
  end
else
  errStr = 'Unknown syntax!';
  irf.log('critical',errStr);	error(errStr);
end

%% Define parameters for selected database
switch omniDatabase
  case 'omni_hour'
    dataSource  = 'omni2';
    dateFormat  = 'utc_yyyymmdd';
    dtMin       = 24*3600;
  case 'omni_min'
    dataSource  = 'omni_min';
    dateFormat  = 'utc_yyyymmddHH';
    dtMin       = 3600;
    dtAtIntervalEnds = 0.5*3600;
  otherwise
    errStr = ['no such source: '  omniDatabase];
    irf.log('critical',errStr);
    error(errStr);
end

%% Define request url and time interval
if isa(tint,'GenericTimeArray') && length(tint)==2
  tintTemp = tint.epochUnix;
  tint = tintTemp(:)';
end

httpRequest = ['https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi?activity=retrieve&spacecraft=' dataSource '&'];
startDate   = irf_time(tint(1)        ,dateFormat);
endDate     = irf_time(tint(2) + dtMin,dateFormat);

%% Define variables to be requested
ii     = strfind(parameter,',');
iEnd   = [ii-1 length(parameter)];
iStart = [1 ii+1];
vars='';
nVar=0;
for jj=1:length(iStart)
  variable=parameter(iStart(jj):iEnd(jj));
  switch lower(variable)
    case 'b',      varOmni2=8 ;varOmni1min=13;
    case 'avgb',   varOmni2=9 ;varOmni1min=-1;
    case 'blat',   varOmni2=10;varOmni1min=-1;
    case 'blong',  varOmni2=11;varOmni1min=-1;
    case 'bx',     varOmni2=12;varOmni1min=14;
    case 'bxgse',  varOmni2=12;varOmni1min=14;
    case 'bxgsm',  varOmni2=12;varOmni1min=14;
    case 'by',     varOmni2=13;varOmni1min=15;
    case 'bygse',  varOmni2=13;varOmni1min=15;
    case 'bz',     varOmni2=14;varOmni1min=16;
    case 'bzgse',  varOmni2=14;varOmni1min=16;
    case 'bygsm',  varOmni2=15;varOmni1min=17;
    case 'bzgsm',  varOmni2=16;varOmni1min=18;
    case 't',      varOmni2=22;varOmni1min=26;
    case 'n',      varOmni2=23;varOmni1min=25;
    case 'nanp',   varOmni2=27;varOmni1min=-1;
    case 'v',      varOmni2=24;varOmni1min=21;
    case 'vx',     varOmni2=-1;varOmni1min=22;
    case 'vy',     varOmni2=-1;varOmni1min=23;
    case 'vz',     varOmni2=-1;varOmni1min=24;
    case 'vlon',   varOmni2=25;varOmni1min=-1;
    case 'vlat',   varOmni2=26;varOmni1min=-1;
    case 'p',      varOmni2=28;varOmni1min=27;
    case 'e',      varOmni2=35;varOmni1min=28;
    case 'beta',   varOmni2=36;varOmni1min=29;
    case 'ma',     varOmni2=37;varOmni1min=30;
    case 'bsnx',   varOmni2=-1;varOmni1min=34;
    case 'bsny',   varOmni2=-1;varOmni1min=35;
    case 'bsnz',   varOmni2=-1;varOmni1min=36;
    case 'ms',     varOmni2=56;varOmni1min=45;
    case 'ssn',    varOmni2=39;varOmni1min=-1;
    case 'dst',    varOmni2=40;varOmni1min=-1;
    case 'ae',     varOmni2=41;varOmni1min=37;
    case 'al',     varOmni2=52;varOmni1min=38;
    case 'au',     varOmni2=53;varOmni1min=39;
    case 'kp',     varOmni2=38;varOmni1min=-1;
    case 'pc',     varOmni2=51;varOmni1min=44;
    case 'f10.7',  varOmni2=50;varOmni1min=-1;
    case 'imfid',  varOmni2=-1;varOmni1min= 4;
    case 'swid',   varOmni2=-1;varOmni1min= 5;
    case 'ts',     varOmni2=-1;varOmni1min=9;
    case 'rmsts',  varOmni2=-1;varOmni1min=10;
    otherwise,     varOmni2=0 ;varOmni1min=-1;
  end
  if strcmp(dataSource,'omni2')
    if varOmni2>0
      vars=[vars '&vars=' num2str(varOmni2)]; %#ok<AGROW>
      nVar=nVar+1;
    end
  else % datasource omni_min
    if varOmni1min>0
      vars=[vars '&vars=' num2str(varOmni1min)]; %#ok<AGROW>
      nVar=nVar+1;
    end
  end
end

%% Request data
url = [httpRequest 'start_date=' startDate '&end_date=' endDate vars];
disp(['url: ' url]);
try
  if verLessThan('matlab','9.1') % < Version less than R2016b
    % Old Matlab unable to disable certificate
    if isunix
      % Try external program, "curl" or "wget".
      prog = ''; args='';
      reqSoftware = {'wget', 'curl'};
      for ii = 1:length(reqSoftware)
        [status, ~] = system(['command -v ', reqSoftware{ii}, ' >/dev/null 2>&1 || { exit 100; }']);
        if(status == 0)
          prog = reqSoftware{ii};
          break
        end
      end
      if isempty(prog)
        % Neither wget or curl is installed.
        errStr = ['You appear to be running a too old version of Matlab, ', ...
          'and did not locate system program to download files. Unable to ', ...
          'automatically access the HTTPS url: ', url];
        irf.log('critical', errStr);
        error(['Unable to access the HTTPS-only server with OMNI data, ',...
          'please consider upgrading Matlab or downloading data manually.']);
      else
        % Download using wget or curl
        % Replace "&" with "\&" to avoid expantion problem in various shell
        % environments.
        urlExternal = strrep(url, '&', '\&');
        if strcmp(prog,'wget')
          % Extra arguments to wget (do not check certificate, and output in
          % stdout)
          args = '--no-check-certificate -qO- ';
        elseif strcmp(prog, 'curl')
          % Extra argument to curl (silent progress bar, and output in
          % stdout)
          args = '--insecure -s ';
        end
        [status, c] = system([prog, ' ', args, urlExternal]);
        if status
          % Failed to run "prog" to download url.
          errStr = ['You appear to be running a too old version of Matlab, and program ', ...
            prog, ' failed. Unable to automatically access the HTTPS url: ', url];
          irf.log('critical', errStr);
          error(['Unable to access the HTTPS-only server with OMNI data, ',...
            'please consider upgrading Matlab or downloading data manually.']);
        end
        getDataSuccess = true;
      end
    else
      % Windows system, no wget / curl per default.
      % However PowerShell should probably work on most modern Windows
      % systems, on some versions this may fail as the powershell may have
      % changed its built in functionallity over the years...
      % This has only been tested on a Win 7 machine as of 2017/10/25.
      urlExternal = strrep(url, '&', '^&'); % Escape ampersand with "^".
      cmd = ['powershell -inputformat none $source = ''', ...
        urlExternal, '''; $dest = [System.IO.Path]::GetTempFileName(); ', ...
        '$wc = New-Object System.Net.WebClient; ', ...
        '$wc.DownloadFile($source, $dest); cat $dest; Remove-Item $dest'];
      [status, c] = system(cmd);
      if status
        errStr = ['You appear to be running a too old version of Matlab. ', ...
          'Unable to automatically access the HTTPS url: ', url];
        irf.log('critical', errStr);
        error(['Unable to access the HTTPS-only server with OMNI data, ',...
          'please consider upgrading Matlab or downloading data manually.']);
      end
      getDataSuccess = true;
    end
  else
    % Download from HTTPS using Matlab's webread.
    if verLessThan('matlab','9.2') % Version less than R2017a
      % Set root certificate pem file to empty disables verification, as
      % Matlab versions before R2017a does not include root certificate
      % used by "Let's encrypt".
      webOpt = weboptions('CertificateFilename', '');
    else
      webOpt = weboptions();
    end
    webOpt.Timeout = 20;
    c = webread(url, webOpt);
    getDataSuccess = true;
  end
catch
  getDataSuccess = false;
end


%% Analyze returned data
if getDataSuccess % success in downloading from internet
  cstart=strfind(c,'YEAR'); % returned by omni2 databse
  if isempty(cstart)
    cstart=strfind(c,'YYYY'); % returned by omni_min database
  end
  if isempty(cstart) % no data returned
    irf.log('warning','Can not get OMNI data from internet!');
    display(c)
    f=[];
    return
  end
  cend=strfind(c,'</pre>')-1;
  if strcmp(dataSource,'omni2')
    fmt=['%f %f %f' repmat(' %f',1,nVar)];
    cc=textscan(c(cstart:cend),fmt,'headerlines',1);
    xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
    f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600;
    for jj=1:nVar
      f(:,jj+1)=cc{jj+3};
    end
  else
    fmt=['%f %f %f %f' repmat(' %f',1,nVar)];
    cc=textscan(c(cstart:cend),fmt,'headerlines',1);
    xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
    f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600+cc{4}*60;
    for jj=1:nVar
      f(:,jj+1)=cc{jj+4};
    end
  end
  f = irf_tlim(f,tint); % leave only data within the time interval
  % Remove FILLVAL, however not good solution. TODO: improve, so that
  % only FILLVAL for the corresponding variable is used
  f(f==9999999.)=NaN;
  f(f==999999.9)=NaN;
  f(f==99999.9)=NaN;
  f(f==9999.99)=NaN;
  f(f==999.99)=NaN;
  f(f==999.9)=NaN;
  f(f==999)=NaN;
  f(f==99.99)=NaN;
else % no success in getting data from internet
  irf.log('warning','Can not get OMNI data from internet!');
  f=[];
end
