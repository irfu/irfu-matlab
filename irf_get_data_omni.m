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
%               'E'     - electric field (mV/m)
%               'P'     - flow pressure (nPa)
%               'beta'  - plasma beta
%               'Ma'    - Alfven Mach number
%               'ssn'   - daily sunspot number
%               'dst'   - DST index
%               'f10.7' - F10.7 flux
%               'pc'    - PC index
%               'ae'    - AE index
%               'al'    - AL index
%               'au'    - AL index
% 
% f=IRF_GET_DATA_OMNI(tint,parameter,database) download from specified database 
%
% database:  'omni2'    - 1h resolution OMNI2 data (default)
%            'omni_min' - 1min resolution OMNI data
%
% Examples:
%   tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 12 31 23 59 0])];
%   ff= irf_get_data_omni(tint,'b,bx,bygsm');
%   ff= irf_get_data_omni(tint,'f10.7');
%   ff= irf_get_data_omni(tint,'b','omni_min');

% http://omniweb.gsfc.nasa.gov/html/ow_data.html

% data description of high res 1min,5min omni data set
% http://omniweb.gsfc.nasa.gov/html/omni_min_data.html#4b
% %
% %
% Year			        I4	      1995 ... 2006
% Day			        I4	1 ... 365 or 366
% Hour			        I3	0 ... 23
% Minute			        I3	0 ... 59 at start of average
% ID for IMF spacecraft	        I3	See  footnote D below
% ID for SW Plasma spacecraft	I3	See  footnote D below
% # of points in IMF averages	I4
% # of points in Plasma averages	I4
% Percent interp		        I4	See  footnote A above
% Timeshift, sec		        I7
% RMS, Timeshift		        I7
% RMS, Phase front normal	        F6.2	See Footnotes E, F below
% Time btwn observations, sec	I7	DBOT1, See  footnote C above
% Field magnitude average, nT	F8.2
% Bx, nT (GSE, GSM)		F8.2
% By, nT (GSE)		        F8.2
% Bz, nT (GSE)		        F8.2
% By, nT (GSM)	                F8.2	Determined from post-shift GSE components
% Bz, nT (GSM)	                F8.2	Determined from post-shift GSE components
% RMS SD B scalar, nT	        F8.2	
% RMS SD field vector, nT	        F8.2	See  footnote E below
% Flow speed, km/s		F8.1
% Vx Velocity, km/s, GSE	        F8.1
% Vy Velocity, km/s, GSE	        F8.1
% Vz Velocity, km/s, GSE	        F8.1
% Proton Density, n/cc		F7.2
% Temperature, K		        F9.0
% Flow pressure, nPa		F6.2	See  footnote G below		
% Electric field, mV/m		F7.2	See  footnote G below
% Plasma beta		        F7.2	See  footnote G below
% Alfven mach number		F6.1	See  footnote G below
% X(s/c), GSE, Re		        F8.2
% Y(s/c), GSE, Re		        F8.2
% Z(s/c), GSE, Re		        F8.2
% BSN location, Xgse, Re	        F8.2	BSN = bow shock nose
% BSN location, Ygse, Re	        F8.2
% BSN location, Zgse, Re 	        F8.2
% 
% AE-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
% AL-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
% AU-index, nT                    I6      See World Data Center for Geomagnetism, Kyoto
% SYM/D index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% SYM/H index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% ASY/D index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% ASY/H index, nT                 I6      See World Data Center for Geomagnetism, Kyoto
% PC(N) index,                    F7.2    See World Data Center for Geomagnetism, Copenhagen

if nargin < 3, % database not specified defaulting to omni2
  datasource='omni2';
  dateformat='yyyymmdd';
elseif nargin == 3, % database specified
  if strcmpi(database,'omni2')
    datasource='omni2';
    dateformat='yyyymmdd';
  elseif strcmpi(database,'omni_min') || strcmpi(database,'min')
    datasource='omni_min';
    dateformat='yyyymmddhh';
  else
    irf_log('fcal','Unknown database, using omni2.');
    datasource='omni2';
    dateformat='yyyymmdd';
  end
end
httpreq=['http://omniweb.gsfc.nasa.gov/cgi/nx1.cgi?activity=retrieve&spacecraft=' datasource '&'];
start_date=irf_time(tint(1),dateformat);
end_date=irf_time(tint(2),dateformat);

i=strfind(parameter,',');
iend=[i-1 length(parameter)];
istart=[1 i+1];
vars='';number_var=0;
for jj=1:length(istart)
    variable=parameter(istart(jj):iend(jj));
    switch lower(variable)
        case 'b', var_number_omni2=8;var_number_omni1min=13;
        case 'avgb', var_number_omni2=9;var_number_omni1min=-1;
        case 'blat', var_number_omni2=10;var_number_omni1min=-1;
        case 'blong', var_number_omni2=11;var_number_omni1min=-1;
        case {'bx','bxgse'}, var_number_omni2=12;var_number_omni1min=14;
        case {'by','bygse'}, var_number_omni2=13;var_number_omni1min=15;
        case {'bz','bzgse'}, var_number_omni2=14;var_number_omni1min=16;
        case 'bygsm', var_number_omni2=14;var_number_omni1min=17;
        case 'bzgsm', var_number_omni2=15;var_number_omni1min=18;
        case 't', var_number_omni2=22;var_number_omni1min=26;
        case 'n', var_number_omni2=23;var_number_omni1min=25;
        case 'nanp', var_number_omni2=27;var_number_omni1min=-1;
        case 'v', var_number_omni2=24;var_number_omni1min=21;
        case 'p', var_number_omni2=28;var_number_omni1min=27;
        case 'e', var_number_omni2=35;var_number_omni1min=28;
        case 'beta', var_number_omni2=36;var_number_omni1min=29;
        case 'ma', var_number_omni2=37;var_number_omni1min=30;
        case 'ssn', var_number_omni2=39;var_number_omni1min=-1;
        case 'dst', var_number_omni2=40;var_number_omni1min=-1;
        case 'ae', var_number_omni2=41;var_number_omni1min=37;
        case 'al', var_number_omni2=52;var_number_omni1min=38;
        case 'au', var_number_omni2=53;var_number_omni1min=39;
        case 'kp', var_number_omni2=38;var_number_omni1min=-1;
        case 'pc', var_number_omni2=51;var_number_omni1min=44;
        case 'f10.7', var_number_omni2=50;var_number_omni1min=-1;
        otherwise, var_number_omni2=0;var_number_omni1min=-1;
    end
    if strcmp(datasource,'omni2'),
      if var_number_omni2>0,
        vars=[vars '&vars=' num2str(var_number_omni2)];
        number_var=number_var+1;
      end
    else % datasource omni_min
      if var_number_omni1min>0,
        vars=[vars '&vars=' num2str(var_number_omni1min)];
        number_var=number_var+1;
      end
    end
end

url=[httpreq 'start_date=' start_date '&end_date=' end_date vars];
disp(['url:' url]);
[c,status]=urlread(url);

if status==1, % success in downloading from internet
    cstart=strfind(c,'YEAR'); % returned by omni2 databse
    if isempty(cstart), 
      cstart=strfind(c,'YYYY'); % returned by omni_min database
    end
    if isempty(cstart), % no data returned
      irf_log('fcal','Can not get OMNI data from internet!');
      f=[];
      return
    end
    cend=strfind(c,'</pre>')-1;
    if strcmp(datasource,'omni2')
      fmt=['%f %f %f' repmat(' %f',1,number_var)];
      cc=textscan(c(cstart:cend),fmt,'headerlines',1);
      xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
      f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600;
      for jj=1:number_var,
        f(:,jj+1)=cc{jj+3};
      end
    else
      fmt=['%f %f %f %f' repmat(' %f',1,number_var)];
      cc=textscan(c(cstart:cend),fmt,'headerlines',1);
      xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
      f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600+cc{4}*60;
      for jj=1:number_var,
        f(:,jj+1)=cc{jj+4};
      end
    end
    f(f==9999999.)=NaN;
    f(f==999999.9)=NaN;
    f(f==99999.9)=NaN;
    f(f==9999.99)=NaN;
    f(f==999.99)=NaN;
    f(f==999.9)=NaN;
    f(f==999)=NaN;
    f(f==99.99)=NaN;
else % no success in getting data from internet
    irf_log('fcal','Can not get OMNI data form internet!');
    f=[];
end

