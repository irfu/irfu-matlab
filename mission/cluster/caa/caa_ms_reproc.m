function caa_ms_reproc(fname)
%CAA_MS_REPROC  Reprocess magnetospheric data
%
% caa_ms_reproc(fname)
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
dirs = textread(fname,'%s');

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
  curr_d = dirs{d};
  cd( [BASE_DIR '/' curr_d])
  
  if ~exist('./.caa_sh_interval','file') && ~exist('./.caa_ms_interval','file')
    cl_id = str2double(curr_d(21));
    if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end
    
    irf_log('proc',[ '-- GETTING -- : ' curr_d]);
    if exist('./mERC.mat','file') || exist('./mEDSI.mat','file')
      !rm -f mERC.mat mEDSI.mat mEDSIf.mat
    end
    
    getData(ClusterProc(pwd),cl_id,'dies');
    getData(ClusterProc(pwd),cl_id,'diespec');
    
    [iso_t,dt] = caa_read_interval;
    st = iso2epoch(iso_t);
    getData(ClusterDB,st,dt,cl_id,'v');
    getData(ClusterDB,st,dt,cl_id,'bfgm');
    getData(ClusterDB,st,dt,cl_id,'b');
    getData(ClusterDB,st,dt,cl_id,'edi');
    getData(ClusterProc(pwd),cl_id,'brs');
    getData(ClusterProc(pwd),cl_id,'edi');
    getData(ClusterProc(pwd),cl_id,'wake');
    
    % Create .caa_ms_interval
    fid = fopen('.caa_ms_interval','w');
    if fid<0
      irf_log('save','problem creating .caa_ms_interval')
      cd(old_pwd),return
    end
    count = fprintf(fid,'%s',epoch2iso(date2epoch(now))); fclose(fid);
    if count<=0
      irf_log('save','problem writing to .caa_ms_interval')
      cd(old_pwd), return
    end
    
  else
    irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
  end
end
cd(old_pwd)
