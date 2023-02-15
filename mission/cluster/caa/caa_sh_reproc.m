function caa_sh_reproc(fname)
%CAA_SH_REPROC  Reprocess magnetosheath/sw data
%
% caa_sh_reproc(fname)
%

% Copyright 2007 Yuri Khotyaintsev

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
dirs = textread(fname,'%s'); MP=[];

if exist('./mPlan.mat','file'), load ./mPlan.mat
else, error('No MPlan.mat found')
end

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
  curr_d = dirs{d};
  cd( [BASE_DIR '/' curr_d])
  
  if ~exist('./.caa_sh_interval','file')
    cl_id = str2double(curr_d(21));
    if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end
    
    yyyy = str2double(curr_d(1:4));
    v_s = sprintf('MPauseY%d',yyyy);
    if ~exist(v_s,'var'), error(['Cannot load ' v_s]), end
    eval([ 'MP=' v_s ';'])
    
    [iso_t,dt] = caa_read_interval;
    st = iso2epoch(iso_t); et = st +dt;
    
    if ~isempty( find( MP(:,1)>=st & MP(:,1)<et ,1) ) || ...
        ~isempty( find( MP(:,2)>st & MP(:,2)<=et ,1) ) || ...
        ~isempty( find( MP(:,1)<=st & MP(:,2)>=et ,1) )
      
      
      irf_log('proc',[ '-- GETTING -- : ' curr_d]);
      if exist('./mERC.mat','file') || exist('./mEDSI.mat','file')
        !rm -f mERC.mat mEDSI.mat mEDSIf.mat
      end
      getData(ClusterProc(pwd),cl_id,'ec','correct_sw_wake');
      getData(ClusterProc(pwd),cl_id,'rawspec');
      getData(ClusterProc(pwd),cl_id,'dies');
      getData(ClusterProc(pwd),cl_id,'diespec');
      getData(ClusterProc(pwd),cl_id,'die');
      
      % Create .caa_sh_interval
      fid = fopen('.caa_sh_interval','w');
      if fid<0
        irf_log('save','problem creating .caa_sh_interval')
        cd(old_pwd),return
      end
      count = fprintf(fid,'%s',epoch2iso(date2epoch(now))); fclose(fid);
      if count<=0
        irf_log('save','problem writing to .caa_sh_interval')
        cd(old_pwd), return
      end
      
    else, irf_log('proc',[ '-- INSIDE MP -- : ' curr_d]);
    end
    
  else
    irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
  end
end
cd(old_pwd)
