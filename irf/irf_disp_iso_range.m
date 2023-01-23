function s = irf_disp_iso_range(data,fmt)
%IRF_DISP_ISO_RANGE  display time ranges in ISO format
%
% irf_disp_iso_range(data,fmt)
%
% See also: EPOCH2ISO
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


if nargin<2, fmt = 0; end

ndata = size(data,1);
s = cell(1,ndata);
for i=1:ndata
  s{i} = [epoch2iso(data(i,1),fmt) ' -- ' epoch2iso(data(i,2),fmt)];
end
if nargout==0
  for i=1:ndata
    disp(s{i})
  end
end
if ndata==1, s = s{1}; end