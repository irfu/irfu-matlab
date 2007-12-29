function caa_get_edi(fname)
%CAA_GET_EDI  get EDI data
%
% caa_sh_get_edi(fname)
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
	
	if ~exist('./mEDI.mat','file')
	  cl_id = str2double(curr_d(21));
	  if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end

          [iso_t,dt] = caa_read_interval;
          st = iso2epoch(iso_t);

          irf_log('proc',[ '-- GETTING -- : ' curr_d]);
	  outd = getData(ClusterDB,st,dt,cl_id,'edi');
	  if ~isempty(outd)
	    getData(ClusterDB,st,dt,cl_id,'v')
            getData(ClusterDB,st,dt,cl_id,'b')
	    getData(ClusterDB,st,dt,cl_id,'bfgm')
            getData(ClusterProc(pwd),cl_id,'edi')
	  else
	    irf_log('load','-- NO EDI DATA --');
	  end

	else
	  irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
	end
end
cd(old_pwd)
