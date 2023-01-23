function caa_reproc(fname,l0vars,l1vars,cfiles)
%CAA_REPROC  Reprocess CAA data
%
% CAA_REPROC(FLIST_FILE,L0VARS,L1VARS,CFILES)
%
% FLIST_FILE - file containing list of directories containing
%       .caa_ms/sh_interval files in following formats:
%       with full subdirs YYYY/YYYYMMDD_hhmm/Cn/YYYYMMDD_hhmm
%       full path and no subdirs /data/caa/l1/YYYY/YYYYMMDD_hhmm/Cn
%       (produced by find_data_dirs.sh).
%
% L0VARS - list of L0 variables separated by '|'
% L1VARS - list of L1 variables separated by '|'
% CFILES - list of files to be removed prior to procedsiing separated by
%       '|' and without the '.mat' extension
%
% See also CAA_SH_REPROC, CAA_MS_REPROC
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,4)

if nargin<4, cfiles = ''; end
if nargin==2, l1vars = ''; end

if ~isempty(l0vars), l0vars = tokenize(l0vars,'|'); end
if ~isempty(l1vars), l1vars = tokenize(l1vars,'|'); end
if ~isempty(cfiles), cfiles = tokenize(cfiles,'|'); end

if isempty(l0vars) && isempty(l1vars) && isempty(cfiles)
  disp('NOTHING TO REPROCESS')
  return
end

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
write_caa_reproc=0;

dirs = textread(fname,'%s');
if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
  dir_s = dirs{d};
  if strcmp(dir_s(1:12),BASE_DIR)
    dir_s = dir_s(14:end); % Convert to relative path
  end
  
  if length(dir_s) <= 21 % /data/caa/l1/2001/20010204_0900/C1
    sdirs = dir([BASE_DIR '/' dir_s '/2*_*']);
    if isempty(sdirs), continue, end
    sdirs = struct2cell(sdirs);
    sdirs = sdirs(1,:);
  else
    sdirs = {''};
  end
  
  for sdi = 1:length(sdirs)
    if isempty(sdirs{sdi}), curr_d = dir_s;
    else, curr_d = [dir_s '/' sdirs{sdi}];
    end
    
    cd( [BASE_DIR '/' curr_d])
    
    if ~write_caa_reproc || ~exist('./.caa_reproc','file') && ( exist('./.caa_sh_interval','file') ...
        || exist('./.caa_ms_interval','file') )
      cl_id = str2double(curr_d(21));
      if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end
      
      irf_log('proc',[ '-- GETTING -- : ' curr_d]);
      
      % Clean files
      if ~isempty(cfiles)
        for i=1:length(cfiles)
          if exist(['./' cfiles{i} '.mat'],'file')
            [s,w] = unix(['rm ./' cfiles{i} '.mat']);
            if s ~= 0
              error(['cannot remove ' cfiles{i} '.mat : ' w]);
            else
              irf_log('proc',['removed ' cfiles{i} '.mat'])
            end
          end
        end
      end
      
      % L0 vars
      if ~isempty(l0vars)
        [iso_t,dt] = caa_read_interval;
        st = iso2epoch(iso_t);
        for i=1:length(l0vars)
          %disp(['L0 : ' l0vars{i}])
          getData(ClusterDB,st,dt,cl_id,l0vars{i});
        end
      end
      
      % L1 vars
      if ~isempty(l1vars)
        for i=1:length(l1vars)
          %disp(['L1 : ' l1vars{i}])
          getData(ClusterProc(pwd),cl_id,l1vars{i});
        end
      end
      
      % Create .caa_ms/sh_interval
      if exist('./.caa_sh_interval','file'), lf = '.caa_sh_interval';
      else, lf = '.caa_ms_interval';
      end
      
      fid = fopen(lf,'w');
      if fid<0
        irf_log('save',['problem creating ' lf])
        cd(old_pwd),return
      end
      count = fprintf(fid,'%s',epoch2iso(date2epoch(now))); fclose(fid);
      if count<=0
        irf_log('save',['problem writing to ' lf])
        cd(old_pwd), return
      end
      
      if write_caa_reproc
        lf = '.caa_reproc'; %#ok<UNRCH>
        fid = fopen(lf,'w');
        if fid<0
          irf_log('save',['problem creating ' lf])
          cd(old_pwd),return
        end
        count = fprintf(fid,'%s',epoch2iso(date2epoch(now))); fclose(fid);
        if count<=0
          irf_log('save',['problem writing to ' lf])
          cd(old_pwd), return
        end
      end
      
    else
      irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
    end
  end
end
cd(old_pwd)
