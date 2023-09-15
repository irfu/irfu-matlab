function caa_deltaoff_batch(fname)
%CAA_DELTAOFF_BATCH  invert previously applied delta offsets
%
% caa_deltaoff_batch(fname)
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

old_pwd = pwd;
dirs = textread(fname,'%s');

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
  s = dirs{d};
  cd (s);
  cl_id = str2double( s(34) );
  [ok, Del] = c_load('D?p12p34',cl_id);
  if ~ok || ~isreal(Del(1)), continue, end

  irf_log('proc',['inverting delta in ' s(33:end)])

  [ok, p12] = c_load('diEs?p12',cl_id);
  if ~ok, error('cannot load p12'), end

  [ok, p34] = c_load('diEs?p34',cl_id);
  if ~ok, error('cannot load p34'), end

  % apply delta offset to p34 instead
  p12(:,2:3) = p12(:,2:3) +ones(size(p12(:,1),1),1)*Del; %#ok<NASGU>
  p34(:,2:3) = p34(:,2:3) +ones(size(p34(:,1),1),1)*Del; %#ok<NASGU>
  Del = -Del*1i; %#ok<NASGU>

  c_eval('D?p12p34 = Del; diEs?p12 = p12; diEs?p34 = p34; save mEDSI.mat D?p12p34 diEs?p12 diEs?p34 -append',cl_id)
  getData(ClusterProc(pwd),cl_id,'die');
end
cd(old_pwd)
