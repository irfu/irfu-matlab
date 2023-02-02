function caa_rerun_L3E(fname)
%CAA_RERUN_L3E  rerun spinfits
%
% caa_rerun_L3E(fname)
%

% Copyright 2006 Yuri Khotyaintsev

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
cd (BASE_DIR)
dirs = textread(fname,'%s');

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
  curr_d = dirs{d};
  cd( [BASE_DIR '/' curr_d])
  
  if ~exist('./.caa_rerun_L3E','file')
    cl_id = str2double(curr_d(21));
    if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end
    
    irf_log('proc',[ '-- RERUNNING -- : ' curr_d]);
    getData(ClusterProc(pwd),cl_id,'dies')
    
    % Create .caa_rerun_L3E
    fid = fopen('.caa_rerun_L3E','w');
    if fid<0
      irf_log('save','problem creating .caa_rerun_L3E')
      cd(old_pwd),return
    end
    count = fprintf(fid,'%s',epoch2iso(date2epoch(now))); fclose(fid);
    if count<=0
      irf_log('save','problem writing to .caa_rerun_L3E')
      cd(old_pwd), return
    end
  else
    irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
  end
end
cd(old_pwd)
