function res = c_load(vs,sp,st,dt,db,dp)
%function c_load(vs,sp,st,dt,db,dp)
%function c_load(vs,sp,st,dt,cdb)
% return O if OK.

error(nargchk(1,5,nargin))

d = c_desc(vs);

if nargin < 2, sp = '.';, end
old_pwd = pwd;
cd(sp);

%try to load from file
if exist([d.file '.mat'],'file'), eval(['load -mat ' d.file ' ' vs]), end

%try to load with getData
if ~exist(vs,'var')
	if strcmp(d.file,'mEdB')
		getData(ClusterProc,d.cl_id,d.quant);
	else
		if (isa(db,'ClusterDB')) 
			cdb = db; 
			cdb = set(cdb,'sp',sp);
		else
			if nargin < 6, dp = '/data/cluster', end
			cdb = ClusterDB(db,dp,sp);
		end
		getData(cdb,st,dt,d.cl_id,d.quant);
	end
	%try again to load from file
	if exist([d.file '.mat'],'file'), c_eval(['load -mat ' d.file ' ' vs]), end
end

%return the result
if exist(vs,'var')
	assignin('caller',vs,eval(vs));
	if nargout>0, res = 0, end
else
	if nargout>0, res = 1, 
	else, c_log('load',['cannot load ' vs])
	end
end

cd(old_pwd)
