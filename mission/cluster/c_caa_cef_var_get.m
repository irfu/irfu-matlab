function out = c_caa_cef_var_get(varName,fileName,varargin)
%C_CAA_CEF_GET_DATA read variable from CEF file
%
% out = C_CAA_CEF_VAR_GET(varName,fileName) reads variable varName from cef
% file fileName. The output out is matrix of data type VariableMat.
%
% out = C_CAA_CEF_VAR_GET({varName1,varName2,..},fileName) reads variables
% varName1, varName2, .. from cef file fileName. The output out is cell
% array where each cell element is corresponding variable as data type
% VariableMat.
%
% out = C_CAA_CEF_VAR_GET(cdfVarName,'tint',tint,'stream') streams variable
% from CAA as data type VariableMat.
%
% Example:
%   R = C_CAA_CEF_GET_DATA('sc_r_xyz_gse','CL_SP_AUX__20101231_010001_20101231_010201.cef.gz')
%
% See also: CEF_INIT, CEF_READ, IRF.DATATYPES

%% Check CEFlib status
persistent okCeflib 

if isempty(okCeflib) 
	% check whether ceflib is properly installed
	okCeflib = datastore('irfu_matlab','okCeflib');
end

if ~okCeflib
	irf.log('critical','ceflib is not installed properly!');
	if nargout > 0, out = []; end
	return
end
%% Define defaults

% returnOutputAsCellArray  - set later in code, if true output should be cell array

%% Check inputs

if nargin == 0, help c_caa_cef_var_get; return; end
if nargin == 1
	errStr = 'Should be at least 2 input parameters';
	irf.log('critical',errStr);
	error(errStr);
end
if ischar(fileName)
	if strcmpi(fileName,'tint') % tint specified, get stream from CAA
		assert(nargin == 4 && ischar(varargin{2}) && strcmpi(varargin{2},'stream'),...
			'Syntax not correct');
		if ischar(varargin{1}) % tint in iso format
			tint = irf_time(varargin{1},'utc>tint');
		elseif isnumeric(varargin{1}) ...
				&& numel(varargin{1}) == 2 % tint as vector [tstart tend]
			tint = varargin{1};
		else
			errStr = 'tint not specified in right format';
			irf.log('critical',errStr);
			error(errStr);
		end
		% stream the the file with data
		currentDir = pwd;
		tempDir = tempname;
		mkdir(tempDir);
		cd(tempDir);
		[datasetName,varName]=caa_get_dataset_name(varName);
		caa_download(tint,datasetName{1},'stream'); % TODO: assumes all variables from the same dataset
		cd(['CAA/' datasetName{1}]);
		d=dir('*.cef.gz');
    if isempty(d), out = []; return, end % No data
		cefFile = d.name;
		cef_init();
		cef_read(cefFile);
		c1 = onCleanup(@() cef_close());
		c2 = onCleanup(@() cd(currentDir));
		c3 = onCleanup(@() rmdir(tempDir,'s'));
	else % read file
		cef_init();
		cef_read(fileName);
		c = onCleanup(@() cef_close());
	end
else
	errStr = '2nd input parameter should be string';
	irf.log('critical',errStr);
	error(errStr);
end

% get time
tt=cef_var('time_tags');
tt=tt';
tt=irf_time( cef_date(tt),'datenum>epoch');

% make variable cell array if it is string
if ischar(varName)
	varName = {varName};
	returnOutputAsCellArray = false;
elseif iscellstr(varName)
	returnOutputAsCellArray = true;
else
	irf.log('critical','varName incorrect format');
	error('varName incorrect format');
end
if iscell(varName)
	out = cell(size(varName));
	for iVar = 1:numel(varName)
		temp = cef_var(varName{iVar});
		out{iVar} = [tt double(temp')];
	end
else
	errStr = 'Input varName has unknown type';
	irf.log('critical',errStr);
	error(errStr);
end


% define output
if numel(out) == 1 && ~returnOutputAsCellArray
	out = out{1}; % return only matrix if one variable requested
end

