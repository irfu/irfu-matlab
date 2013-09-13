function f = irf_get_data( tint, parameter , database, format)
%IRF_GET_DATA get data form different databases
%
% f=IRF_GET_DATA(tint,parameter,database) download parameters for specified time interval
%       default f is column vector, first column time and parameters in next columns
%       parameters - string, data names separated by comma (case does not matter).
%
% f=IRF_GET_DATA(parameter,database) checking if tint defined in calling environment and using that one 
%       if tint not defined reading all data
%
% f=IRF_GET_DATA(parameter,database,format) get data in specified format (specific for databases)
%
% database:  'omni2'    - 1h resolution OMNI2 data (default), see also IRF_GET_DATA_OMNI
%            'omni_min' - 1min resolution OMNI data
%            'caa'      - Cluster Active Archive, see also C_CAA_VAR_GET
%                         formats: 'caa','mat','units'
%
% Examples:
%   tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 01 02 23 59 0])];
%   ff= irf_get_data(tint,'b,bx,bygsm','omni_min');
%   ff= irf_get_data(tint,'f10.7','omni');

% http://omniweb.gsfc.nasa.gov/html/ow_data.html
nargs=nargin; % number of defined input arguments
timeIntervalNotDefined = false;
if ischar(tint), % tint not specified
  timeIntervalNotDefined=true;
  if nargin==2,
    database=parameter;parameter=tint;
  elseif nargin==3,
    format=database;database=parameter;parameter=tint;
  end
  irf_log('fcal',['Reading ' parameter ' from database: ' database '.']);
  irf_log('fcal','tint not defined, reading all data.');
%   Possibility to be smart guessing tint (maybe not good idea, needs special flag?)   
%	irf_log('fcal','Time interval not specified. Analyzing if tint exists.');
%   if evalin('caller','exist(''tint'',''var'')'),
%     irf_log('fcal','Using existing tint variable values.');
%     timeIntervalNotDefined=false;
%     tint=evalin('caller','tint');
%     nargs=nargs+1;
%   else
%       irf_log('fcal','tint not defined, reading all data.');
%   end
end
if nargs == 0,
  help irf_get_data;
  return;
end
if ~exist('database','var'),
  irf_log('fcal','Database not defined.');
  return;
end
if ~exist('format','var'),	% if format is not defined
	 format=[];				% then default format is empty
end

if nargout==1, f=[];end % default return empty

switch lower(database)
  case 'omni'
    irf_log('fcal','WARNING!!! Returning 1h OMNI2 data (use database omni_min if needed high res data)')
    f=irf_get_data_omni(tint,parameter,'omni2');
  case 'omni2'
    f=irf_get_data_omni(tint,parameter,database);
  case 'omni_min'
    f=irf_get_data_omni(tint,parameter,database);
  case 'caa'
    if isempty(format),
      format='mat'; % default value
    end
    if timeIntervalNotDefined,
      f=c_caa_var_get(parameter,format);
    else
      f=c_caa_var_get(parameter,format,'tint',tint);
    end
  otherwise
      error(['Unknown database: ' database]);
end
