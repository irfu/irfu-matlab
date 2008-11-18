function irf_disp_iso_range(data,fmt)
%IRF_DISP_ISO_RANGE  display time ranges in ISO format
%
% irf_disp_iso_range(data,fmt)
%
% See also: EPOCH2ISO
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


if nargin<2, fmt = 0; end

for i=1:size(data,1)
	disp([epoch2iso(data(i,1),fmt) ' -- ' epoch2iso(data(i,2),fmt)])
end
