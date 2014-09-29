function idx = irf_latched_idx(sig)
%IRF_CHK_LATCHED_SIG  find ideces of latched data
%
% idx = irf_chk_latched_sig(sig)
%
% Find indices of latched data (3 consecuetive constant values)

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if min(size(sig))>1
  errS = 'wrong input dimension: expecting a vector';
  irf.log('critical',errS), error(errS)
end
if numel(sig)<3, idx = false(size(sig)); return, end
% Data (MMS) will most likely be a row vector, if column vector, flip it.
if size(sig,1)<size(sig,2), sig=sig'; end;

idx1 = diff(sig)==0;
idx2 = (diff(idx1) == 0) & idx1(1:end-1);
newIdx1 = [idx2; false] | [idx2(1); idx2];
idx = [newIdx1; false] | [newIdx1(1); newIdx1];

