function idx = irf_latched_idx(sig, minPoints)
% IRF_LATCHED_IDX  find ideces of latched data
%
% idx = irf_chk_latched_sig(sig, minPoints)
%
% Find indices of latched data ('minPoints' consecuetive constant values)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,2);
if(nargin==1), minPoints=3; end
if min(size(sig))>1
  errS = 'wrong input dimension: expecting a vector';
  irf.log('critical',errS), error(errS)
end
if numel(sig)<minPoints, idx = false(size(sig)); return, end
% Data (MMS) will most likely be a row vector, if column vector, flip it.
if size(sig,1)<size(sig,2), sig=sig'; end

% idx1 = diff(sig)==0;
% idx2 = (diff(idx1) == 0) & idx1(1:end-1);
% newIdx1 = [idx2; false] | [idx2(1); idx2];
% idx = [newIdx1; false] | [newIdx1(1); newIdx1];

idx = false(size(sig)); % Pre allocate to false (ie. not stuck).
% New code, look for at least minPoints stuck points.
if(minPoints>=10) % FIXME, This limit may need updating..
  % This algorithm is quick when looking for large minPoints but becomes
  % rather slow as the number of segments (with data stuck>=minPoints)
  % increases. For MMS burst mode, locating times when probe is stuck for
  % say one second would require minPoints to be 8192 which is rather slow
  % in the second algorithm below.
  ind1 = [0; find(diff(sig)~=0)]+1;
  % Length of each segment
  lengthSeg = [diff(ind1); length(sig)-ind1(end)+1];
  % Any segments stuck for more than minPoints
  ind2 = find(lengthSeg>=minPoints);
  % For each segment stuck for long enough, change idx to true.
  for ii=1:length(ind2)
    idx(ind1(ind2(ii)):ind1(ind2(ii))+lengthSeg(ind2(ii))-1) = true;
  end
else
  % This algorithm is quick when looking for small minPoints but becomes
  % rather slow as minPoints increases.
  % Loop through and ensure diff(sig)==0 for at least minPoints, then build
  % up the index again. This will be slow when minPoints is large, and
  % quick when minPoints is small, regardless of number of segments found.
  % Break down
  for ii=1:minPoints-1
    if(ii==1)
      ind1 = diff(sig)==0;
    else
      ind1 = diff(ind1)==0 & ind1(1:end-1);
    end
  end
  % Build up
  ind2 = [ind1; false] | [ind1(1); ind1];
  for ii=1:minPoints-2
    ind2 = [ind2; false] | [ind2(1); ind2];
  end
  idx = ind2;
end

end