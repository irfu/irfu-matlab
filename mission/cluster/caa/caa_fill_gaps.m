function [data, filldata] = caa_fill_gaps(data,te,max_gaps)
%CAA_FILL_GAPS(data,te)  fill gaps in the of a dataset
%
% res = caa_fill_gaps(data,te)
%       append NaNs to the end of DATA untill TE
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,3)

if size(data,1)<2
  irf_log('proc','cannot fill gaps (not enough samples)')
  return
end

filldata = [];

fs = c_efw_fsample(data);
if fs<=0
  irf_log('proc','cannot fill gaps')
  return
end

ngap = fix( (te -data(end,1))*fs - 1); % NOTE: Possible error here. Not filling the last point up to TE,
%       if TE-last_point < 2/fs !   (ML)

if nargin > 2
  if ngap > max_gaps
    irf_log('proc',sprintf('Request to fill %d gaps at %s gaps exceeds max_gap.',ngap,epoch2iso(data(end,1),1)))
    return
  end
end

if ngap>0
  tt = zeros(ngap+1,size(data,2));
  tt(:,1) = linspace(data(end,1),data(end,1)+ngap/fs,ngap+1);
  tt = tt(2:end,:);
  tt(:,2:end) = NaN;
  if nargout < 2
    data = [data; tt];
  elseif nargout == 2
    filldata = tt;
  end
  % Do not report 1 point gaps
  if ngap>1
    irf_log('proc',...
      sprintf('filling %d gaps at %s',ngap,epoch2iso(data(end,1),1)))
  end
end
