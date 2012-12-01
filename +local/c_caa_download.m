function out=c_caa_download(varargin)
% LOCAL.C_CAA_DOWNLOAD download files to /data/caa
%
%   LOCAL.C_CAA_DOWNLOAD(dataset) download all dataset
%
% 	See also CAA_DOWNLOAD, IRF_FUNCTION_B.
%

% $Id$

if nargin==1 && ischar(varargin{1})
	dataset=varargin{1};
else
	irf_log('fcal','See syntax: help local.c_caa_download');
end

irf_log('dsrc','Checking list of available times');
tt=caa_download(['list:' dataset]);
if numel(tt)==0, 
	disp('Dataset does not exist or there are no data');
	return;
else
	irf_log('dsrc',['Checking inventory: ' irf_time(tt.TimeInterval(1,:),'tint2iso')]);
	ttInventory = caa_download(tt.TimeInterval(1,:),['list:' dataset]);
end

out=ttInventory;
