function caa_sh_get_cis(fname)
%CAA_SH_GET_CIS  get cis data
%
% caa_sh_get_cis(fname)
%

% Copyright 2007 Yuri Khotyaintsev

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
dirs = textread(fname,'%s');
force_reload=0;

c_ctl('init')
c_ctl('set',5,'isdat_db','130.238.30.32:8')

if exist('/data/caa/l1/mPlan.mat','file'), load /data/caa/l1//mPlan.mat
else, error('No MPlan.mat found')
end

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
  curr_d = dirs{d};
  cd( [BASE_DIR '/' curr_d])

  if force_reload || ~exist('./mCIS.mat','file')
    cl_id = str2double(curr_d(21));
    if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end

    yyyy = str2double(curr_d(1:4));
    v_s = sprintf('MP%dY%d',cl_id,yyyy);
    if ~exist(v_s,'var'), error(['Cannot load ' v_s]), end

    eval([ 'MP=' v_s ';'])

    [iso_t,dt] = caa_read_interval;
    st = iso2epoch(iso_t); et = st +dt;

    if ~isempty( find( MP(:,1)>=st & MP(:,1)<=et ) ) || ...
        ~isempty( find( MP(:,2)>=st & MP(:,2)<=et ) ) || ...
        ~isempty( find( MP(:,1)<=st & MP(:,2)>=et ) )


      irf_log('proc',[ '-- GETTING -- : ' curr_d]);
      getData(ClusterDB,st,dt,cl_id,'sax');
      getData(ClusterDB,st,dt,cl_id,'b');
      getData(ClusterDB,st,dt,cl_id,'vcis');
      getData(ClusterDB,st,dt,cl_id,'v');
      %            getData(ClusterProc(pwd),cl_id,'vce');

    else, irf_log('proc',[ '-- INSIDE MP -- : ' curr_d]);
    end
  else
    irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
  end
end
cd(old_pwd)
