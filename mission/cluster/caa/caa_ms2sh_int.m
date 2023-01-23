function caa_ms2sh_int(sc_list)
%CAA_MS2SH_INT  change MS interval to SH
%
% CAA_MS2SH_INT(SC_LIST)
%     Run the following processing:
%          'ec','correct_sw_wake'
%          'rawspec'
%          'dies'
%          'diespec'
%          'die' (if mEDSIf.mat exists)
%     Then zero PSWAKE and LOWAKE, and finally create .caa_sh_interval
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<1, sc_list = 1:4; end

old_pwd = pwd;

for cl_id=sc_list
  cdir = [old_pwd '/C' num2str(cl_id)];
  if ~exist(cdir, 'dir'), continue, end
  
  d = dir([cdir '/2*_*']);
  if isempty(d), continue, end
  
  for jj=1:length(d)
    curdir = [cdir '/' d(jj).name];
    if ~exist([curdir '/.interval'],'file') || ...
        exist([curdir '/.caa_sh_interval'],'file')
      irf_log('proc',[ '-- SKIPPING -- : C' num2str(cl_id) '/' d(jj).name]);
      continue
    end
    cd(curdir)
    
    irf_log('proc',[ '-- GETTING -- : C' num2str(cl_id) '/' d(jj).name]);
    
    getData(ClusterProc(pwd),cl_id,'ec','correct_sw_wake');
    getData(ClusterProc(pwd),cl_id,'rawspec');
    getData(ClusterProc(pwd),cl_id,'dies');
    getData(ClusterProc(pwd),cl_id,'diespec');
    if exist('./mEDSIf.mat','file')
      getData(ClusterProc(pwd),cl_id,'die');
    end
    
    % Take care of WAKES
    for pp=[12 34]
      for wa={'LO','PS'}
        v_s = sprintf('%sWAKE%dp%d',wa{:},cl_id,pp);
        [ok,ttt] = c_load(v_s); %#ok<NASGU>
        if ok
          eval([v_s '=[]; save mEFW.mat ' v_s ' -APPEND'])
          irf_log('proc',[v_s,' -> []'])
        end
      end
    end
    
    !rm .caa_ms_interval
    
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
    
    cd(old_pwd)
  end
end
