function res = caa_comp_noise_spec_batch(fname,probe_p)

if nargin<2, probe_p = 34; end
count = 0;

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';

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

    cl_id = str2double(curr_d(21));
    if isnan(cl_id) || cl_id>4 || cl_id<1, error(['wrong directory ' curr_d]), end

    irf_log('proc',[ '-- GETTING -- : ' curr_d]);

    [iso_t,dt] = caa_read_interval;
    if isempty(iso_t), continue, end
    st = iso2epoch(iso_t);

    [ok,tmode] = c_load('mTMode?',cl_id);

    if ~ok || ~any(tmode)
      irf_log('proc',[ '-- SKIPPING -- : ' curr_d]);
      continue
    end

    ok = 0; data = [];
    for i=1:length(probe_p)
      if ~ok || isempty(data)
        [ok,data] = c_load(sprintf('wE?p%d', probe_p(i)),cl_id);
      end
    end
    if ~ok || isempty(data), continue, end

    [ok,whip] = c_load('WHIP?',cl_id);

    if ~ok, continue, end

    try
      ttt = caa_comp_noise_spec(data,whip);
      count = count + 1;
      res(count,1) = st;
      res(count,2) = ttt;
    catch
      continue
    end
  end
end

cd (old_pwd)





