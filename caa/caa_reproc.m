function caa_reproc(fname)
%CAA_REPROC  Reprocess L2/L3 E data
%
% caa_reproc(fname)
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
dirs = textread(fname,'%s');

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
	curr_d = dirs{d};
	cd( [BASE_DIR '/' curr_d])
	
	if ~exist('./.caa_repr_interval','file')
		cl_id = str2double(curr_d(21));
		if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end
        
		irf_log('proc',[ '-- GETTING -- : ' curr_d]);
		if exist('./mERC.mat','file') || exist('./mEDSI.mat','file')
			!rm -f mERC.mat mEDSI.mat mEDSIf.mat
		end
		getData(ClusterProc(pwd),cl_id,'dies');
		getData(ClusterProc(pwd),cl_id,'diespec');
		getData(ClusterProc(pwd),cl_id,'die');

		% Create .caa_repr_interval
		fid = fopen('.caa_repr_interval','w');
		if fid<0
			irf_log('save','problem creating .caa_repr_interval')
			cd(old_pwd),return
		end
		count = fprintf(fid,'%s',epoch2iso(date2epoch(now))); fclose(fid);
		if count<=0
			irf_log('save','problem writing to .caa_repr_interval')
			cd(old_pwd), return
		end
			
	else
		irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
	end
end
cd(old_pwd)
