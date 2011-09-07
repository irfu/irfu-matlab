function f = irf_get_data( tint, parameter , database, format)
%IRF_GET_DATA get data form different databases
%
% f=IRF_GET_DATA(tint,parameter,database) download parameters for specified time interval
%       default f is column vector, first column time and parameters in next columns
%       parameters - string, data names separated by comma (case does not matter).
%
% database:  'omni2'    - 1h resolution OMNI2 data (default), see also IRF_GET_DATA_OMNI
%            'omni_min' - 1min resolution OMNI data
%            'caa'      - Cluster Active Archive, see also C_CAA_VAR_GET
% Examples:
%   tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 01 02 23 59 0])];
%   ff= irf_get_data(tint,'b,bx,bygsm','omni_min');
%   ff= irf_get_data(tint,'f10.7','omni');

% http://omniweb.gsfc.nasa.gov/html/ow_data.html

if nargin == 0,
  help irf_get_data;
  return;
elseif nargin < 3;
  irf_cal('fcal','Wrong number of inputs, see help.');
  return;
elseif nargin==3,
  format=[]; % default format is empty
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
      f=c_caa_var_get(parameter,format,'tint',tint);
  otherwise
      error(['Unknown database: ' database]);
end