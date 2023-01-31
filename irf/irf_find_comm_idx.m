function [ii1,ii2]=irf_find_comm_idx(d1,d2,threshold)
%IRF_FIND_COMM_IDX   Find common index of two time series
%
% [IDX1,IDX2] = IRF_FIND_COMM_IDX(E1,E2,[THRESHOLD])
%
% This function can be used when combining two datasets with
% gaps, like p12 and p34. The timeslines are assumed to be identical, if
% the difference is below THERESHOLD (default 1e-5)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin==2, threshold = 1e-5; end

DEBUG = 0;

ii1 = []; ii2 = [];

% we use only time
t1 = d1(:,1); t2 = d2(:,1);
s1 = 1; s2 = 1; % start index
if min(t1(2:end)-t1(1:end-1))<0 % check for negative time jumps
  t1 = rem_neg_time_jumps(t1);
end
if min(t2(2:end)-t2(1:end-1))<0 % check for negative time jumps
  t2 = rem_neg_time_jumps(t2);
end

n_loop = 0;
while 1
  la = min(length(t1(s1:end)),length(t2(s2:end))) - 1;
  d = t1(s1:s1+la)-t2(s2:s2+la);
  ii = find(abs(d) > threshold);
  
  if isempty(ii) % time lines are identical
    % save the interval and break
    ii1 = [ii1 s1:s1+la]; ii2 = [ii2 s2:s2+la];
    if DEBUG, disp('break at 1'), end %#ok<UNRCH>
    break
  end
  
  ii = ii(1); % first point where times differ
  if ii>1
    % save the interval and shift
    if DEBUG, disp('save interval and shift'), end %#ok<UNRCH>
    ii1 = [ii1 s1:s1+ii-2]; ii2 = [ii2 s2:s2+ii-2];
    s1 = s1 + ii - 1; s2 = s2 + ii - 1;
  end
  if d(ii) > 0
    % shift t2
    if DEBUG, disp('shift t2'), end %#ok<UNRCH>
    
    i_t0 = find(abs(t2(s2:end)-t1(s1))<threshold);
    if isempty(i_t0)
      if t1(s1) > t2(end)
        if DEBUG, disp('break at 2'), end %#ok<UNRCH>
        break
      end
      s1 = s1 + 1;
      if s1>length(t1)
        if DEBUG, disp('break at 21'), end %#ok<UNRCH>
        break
      end
      it = find(t2>=t1(s1));
      if isempty(it)
        if DEBUG, disp('break at 22'), end %#ok<UNRCH>
        break
      end
      s2 = it(1); clear it
    else
      s2 = s2 + i_t0(1) - 1;
    end
  else
    %shift t1
    if DEBUG, disp('shift t1'), end %#ok<UNRCH>
    i_t0 = find(abs(t1(s1:end)-t2(s2)) < threshold);
    if isempty(i_t0)
      if t2(s2) > t1(end)
        if DEBUG, disp('break at 3'), end %#ok<UNRCH>
        break
      end
      s2 = s2 + 1;
      if s2>length(t2)
        if DEBUG, disp('break at 31'), end %#ok<UNRCH>
        break
      end
      it = find(t1>=t2(s2));
      if isempty(it)
        if DEBUG, disp('break at 32'), end %#ok<UNRCH>
        break
      end
      s1 = it(1); clear it
    else
      s1 = s1 + i_t0(1) - 1;
    end
  end
  n_loop = n_loop + 1;
  if DEBUG, fprintf('gap # %d\n',n_loop), end %#ok<UNRCH>
end
end

function tv=rem_neg_time_jumps(tv)
curralen=length(tv);
i=1;
while i<curralen
  if tv(i)>tv(i+1)
    tv(i+1)=[];
    curralen=curralen-1;
    i=i-1;
  end
  i=i+1;
end
end
