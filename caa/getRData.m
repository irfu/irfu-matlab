function getRData(sp,st,dt,sc_list,cdb)
%getRData prepare data from all SC: do all steps until spinfits and despin
% getRData(sp,start_time,dt,[sc_list],[cdb])
% Input:
% sp - storage directory
% start_time - ISDAT epoch
% dt - length of time interval in sec
% sc_list - list of SC [optional]
% cdb - CLusterDB object [optional] 
%	// default: ClusterDB('ice:10|disco:10','/data/cluster')
% 
% Example:
% getRData('/home/yuri/caa-data/20020304', ...
%	toepoch([2002 03 04 10 00 00]),30*60)
%
% $Id$
%
% See also ClusterDB/getData, ClusterProc/getData

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(3,5,nargin))

if nargin<4, sc_list=1:4; end
if nargin<5, cdb = ClusterDB('ice:10|disco:10','/data/cluster',sp); end

vars = {'e','p','a','sax','whip','b','edi','ncis','vcis','vce'};

for cl_id=sc_list
	for k=1:length(vars)
		getData(cdb,st,dt,cl_id,vars{k});
	end
end

cp=ClusterProc(sp);

vars = {'dies','die'};
for cl_id=sc_list
	for k=1:length(vars)
		getData(cp,cl_id,vars{k});
	end
end

