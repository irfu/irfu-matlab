function caa_export_batch(cl_id,outdir,m_vars,m_vers)
%CAA_EXPORT_BATCH run CAA_EXPORT in a batch script
%
% caa_export_batch(cl_id,sdir,m_vars,m_vers)
%
% $Id$

% Copyright 2004-2006 Yuri Khotyaintsev

if length(m_vars) ~= length(m_vars)
	error('M_VARS and M_VERS must have the same number of elements')
end

%{ 
% Load Quality information
if exist('./mInfo.mat','file')
	load -mat mInfo caa_q
end
if ~exist('caa_q','var')
	disp('cannot load quality information from mInfo.mat')
	disp('please run_caa_quality')
	return
end
%}

QUALITY = 3;

sp = pwd;
cd(outdir)
for j=1:length(m_vars)
	v = m_vars{j};
	vers = str2num(m_vers{j});
	vers_s = num2str(vers);
	if vers<10, vers_s = ['0' vers_s]; end
	lev = str2num(v(2));
	caa_vs = v(4:end);
	caa_export(lev,caa_vs,cl_id,QUALITY,vers_s,sp)
end
cd(sp)
