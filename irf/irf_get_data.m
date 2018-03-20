function f = irf_get_data(varargin)
%IRF_GET_DATA get data form different databases
%
% f = IRF_GET_DATA(tint,parameter,database) download parameters for
% specified time interval. Default format is f being column vector where
% first column is time and parameters are in next columns. In future it
% will change to TSeries being the default returned format.
%       parameters - string, data names separated by comma (case does not matter).
%
% f=IRF_GET_DATA(parameter,database) checking if tint defined in calling
%   environment and using that one, if tint not defined reading all data.
%
% f=IRF_GET_DATA(parameter,database,format) get data in specified format
%    which can be specific for different databases.
%
% database:  'omni2'    - 1h resolution OMNI2 data (default), see IRF_GET_DATA_OMNI
%            'omni_min' - 1min resolution OMNI data, see IRF_GET_DATA_OMNI
%            'caa' or 'csa' - Cluster Science Archive, see C_CAA_VAR_GET
%
% Examples:
%   tint=[irf_time([2006 01 01 1 1 0]) irf_time([2006 01 02 23 59 0])];
%   ff= irf_get_data(tint,'b,bx,bygsm','omni_min');
%   ff= irf_get_data(tint,'f10.7','omni');

% https://omniweb.gsfc.nasa.gov/html/ow_data.html
nargs=nargin; % number of defined input arguments
if nargs == 0
  help irf_get_data;
  return;
elseif nargs >= 2 && ischar(varargin{1}) % tint not specified
  tint = [];
	parameter=varargin{1};
	database=varargin{2};
	format=varargin(3:end);
  irf.log('notice',['Reading ' parameter ' from database: ' database '.']);
  irf.log('notice','tint not defined, reading all data.');
elseif nargs >= 3 
	tint = varargin{1};
	parameter=varargin{2};
	database=varargin{3};
	format=varargin(4:end);
end

switch lower(database)
  case 'omni'
    irf.log('warning','Returning 1h OMNI2 data (use database omni_min if needed high res data)')
    f=irf_get_data_omni(tint,parameter,'omni2');
  case 'omni2'
    f=irf_get_data_omni(tint,parameter,database);
  case 'omni_min'
    f=irf_get_data_omni(tint,parameter,database);
  case 'caa'
		f=c_caa_var_get(parameter,format{:},'tint',tint);
  otherwise
      error(['Unknown database: ' database]);
end
