function  [onoff,cover] = caa_read_coverage(path_s)
%CAA_READ_COVERAGE  read on/off and coverage files
%
% [onoff,cover] = caa_read_coverage([path_s])
%

% Copyright 2005 Yuri Khotyaintsev

old_pwd = pwd;
if nargin<1, path_s = pwd; end
cd(path_s)

if ~exist('EFWONOFF_COMM.dat','file') || ~exist('COVERAGE_COMM.dat','file')
  error('COVERAGE_COMM.dat and/or EFWONOFF_COMM.dat are not found')
end

[iso_t,cl_id_onoff,on_off]=textread('EFWONOFF_COMM.dat','%s %d %d',-1);
t = iso2epoch(cell2mat(iso_t));

[iso_ts,iso_te,cl_id_cov,tmmode]=textread('COVERAGE_COMM.dat','%s %s %d %d',-1);

ts = iso2epoch(cell2mat(iso_ts));
te = iso2epoch(cell2mat(iso_te));

for cli=1:4
  ii = find(cl_id_onoff==cli);
  onoff(cli) = {[t(ii) on_off(ii)]};
  ii = find(cl_id_cov==cli);
  cover(cli) = {[ts(ii) te(ii) tmmode(ii)]};
end
cd(old_pwd)
