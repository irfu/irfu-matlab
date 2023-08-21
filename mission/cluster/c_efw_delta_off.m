function dof = c_efw_delta_off(data,cl_id)
%C_EFW_DELTA_OFF get delta offsets from database
%
% delta_off = c_efw_delta_off(t,cl_id)
%
% Return precomputed delta (p12 vs p34) offsets
%
% See also CAA_COROF_DELTA, CAA_UPDATE_DELTAOFF
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

DT = 86400*3; % 1/2 of DT used in caa_update_deltaoff

da = [];
load caa/deltaoff.mat
c_eval('da=DeltaOff?;',cl_id);

if (data(end,1) < da(1,1)-DT) || (data(1,1) > da(end,1)-DT) || all(isnan(data(:,1)))
  irf_log('proc','No precomputed Delta offsets for this time')
  dof = [];
  return
end

if length(data(:,1))==1
  t = data(1,1);
  ii1 = find( da(:,1) >= t-DT );
  ii2 = find( da(:,1) < t+DT );
  ii1 = ii1(1); ii2 = ii2(end);
  if any(isnan(da(ii1,2:3))) && ~any(isnan(da(ii2,2:3)))
    dof = da(ii2,2:3);
    irf_log('proc',sprintf('C%d: delta off from the previous interval',cl_id))
  elseif ~any(isnan(da(ii1,2:3))) && any(isnan(da(ii2,2:3)))
    dof = da(ii1,2:3);
    irf_log('proc',sprintf('C%d: delta off from the next interval',cl_id))
  elseif any(isnan(da(ii1,2:3))) && any(isnan(da(ii2,2:3)))
    dof = [];
    irf_log('proc',sprintf('C%d: no Delta offsets',cl_id))
  else
    da = da(~isnan(da(:,2)),:);
    dof = interp1(da(:,1),da(:,2:3),t,'spline');
  end

  % This is a hack. It seems here that p34 is drifting away as p12
  % correlates much better with the other SC.
  if cl_id==2 && t>=toepoch([2003 9 1 0 0 0]) && t<toepoch([2003 12 1 0 0 0])
    dof(2) = -dof(2)*1i;
  end
else
  dof = interp1(da(:,1),da(:,2:3),data(:,1),'spline');
  dof = [data(:,1) dof];
end


