function f = omni( tint, parameter , database)
%OMNI download omni data
%
% DEVELOPMENT VERSION !!! maybe put everything in universal routine irf_get_data
%
% f=OMNI(tint,parameter) download parameters for specified time interval
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
%               'P'     - flow pressure (nPa)
%               'beta'     - plasma beta
%               'ssn'   - daily sunspot number
%               'dst'   - DST index
%               'f10.7' - F10.7 flux
%               'pc'    - PC index
%               'ae'    - AE index
%               'al'    - AL index
%               'au'    - AL index
% 
% f=OMNI(tint,parameter,database) download from specified database 
%
% database:  'omni2'    - 1h resolution OMNI2 data (default)
%            'omni_min' - 1min resolution OMNI data
%
% Examples:
%   tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 12 31 23 59 0])];
%   ff= omni(tint,'b,bx,bygsm');
%   ff= omni(tint,'f10.7');
%   ff= omni(tint,'b','omni_min');

% http://omniweb.gsfc.nasa.gov/html/ow_data.html

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
        case 'b', var_number=8;
        case 'avgb', var_number=9;
        case 'blat', var_number=10;
        case 'blong', var_number=11;
        case {'bx','bxgse'}, var_number=12;
        case {'by','bygse'}, var_number=13;
        case {'bz','bzgse'}, var_number=14;
        case 'bygsm', var_number=14;
        case 'bzgsm', var_number=15;
        case 't', var_number=22;
        case 'n', var_number=23;
        case 'nanp', var_number=27;
        case 'v', var_number=24;
        case 'p', var_number=28;
        case 'beta', var_number=36;
        case 'ssn', var_number=39;
        case 'dst', var_number=40;
        case 'ae', var_number=41;
        case 'al', var_number=52;
        case 'au', var_number=53;
        case 'kp', var_number=38;
        case 'pc', var_number=51;
        case 'f10.7', var_number=50;
        otherwise, var_number=0;
    end
    if var_number>0,
        vars=[vars '&vars=' num2str(var_number)];
        number_var=number_var+1;
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
      irf_log('fcal','Can not get OMNI data form internet!');
      f=[];
      return
    end
    cend=strfind(c,'</pre>')-1;
    fmt=['%f %f %f' repmat(' %f',1,number_var)];
    cc=textscan(c(cstart:cend),fmt,'headerlines',1);
    xx=double([cc{1} repmat(cc{1}.*0+1,1,2) repmat(cc{1}.*0,1,3)]);
    f(:,1)=irf_time(xx)+(cc{2}-1)*3600*24+cc{3}*3600;
    for jj=1:number_var,
        f(:,jj+1)=cc{jj+3};
    end
    f(f==9999.9)=NaN;
else % no success in getting data from internet
    irf_log('fcal','Can not get OMNI data form internet!');
    f=[];
end

