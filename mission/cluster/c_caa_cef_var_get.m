function out = c_caa_cef_var_get(varName,fileName)
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
% Example:
%   R = CEF_GET_DATA('sc_r_xyz_gse','CL_SP_AUX__20101231_010001_20101231_010201.cef.gz')
%
% See also: CEF_INIT, CEF_READ, IRF.DATATYPES

% read file
cef_init();
cef_read(fileName);

% get time 
tt=cef_var('time_tags');
tt=tt';
tt=irf_time( cef_date(tt),'datenum2epoch');

% get variable
if ischar(varName)
	varName = {varName};
end
if iscell(varName)
	out = cell(size(varName));
	for iVar = 1:numel(varName)
		temp = cef_var(varName);
		out{iVar} = [tt double(temp')];
	end
else
	errStr = 'Input varName has unknown type';
	irf.log('critical',errStr);
	error(errStr);
end


% define output
if numel(out) == 1,
	out = out{1}; % return only matrix if one variable requested
end

% close file
cef_close();
end

