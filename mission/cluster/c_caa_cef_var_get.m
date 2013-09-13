function out = c_caa_cef_var_get(varName,fileName)
%CEF_GET_DATA read variable form CEF file
%
% Example:
%   R = CEF_GET_DATA('sc_r_xyz_gse','CL_SP_AUX__20101231_010001_20101231_010201.cef.gz')
%
% See also: CEF_INIT, CEF_READ

% read file
cef_init();
cef_read(fileName);

% get time 
tt=cef_var('time_tags');
tt=tt';
tt=irf_time( cef_date(tt),'datenum2epoch');

% get variable
tempR = cef_var(varName);
tempR = tempR';

% define output
out = [tt double(tempR)];

% close file
cef_close();
end

