function dof = c_efw_delta_off(data,cl_id)
%C_EFW_DELTA_OFF get delta offsets from database
%
% delta_off = c_efw_delta_off(t,cl_id)
%
% Return precomputed delta offsets
%
% See also CAA_UPDATE_DELTAOFF
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

DT = 86400*3; % 1/2 of DT used in caa_update_deltaoff

da = [];
load caa/deltaoff.mat
c_eval('da=DeltaOff?;',cl_id);

if (data(end,1) < da(1,1)-DT) || (data(1,1) > da(end,1)-DT)
	irf_log('proc','No precomputed Delta offsets for this time')
	dof = [];
	return
end

if length(data(:,1))==1
	t = data(1,1);
	ii1 = find( da(:,1) >= t );
	ii2 = find( da(:,1) < t );
	ii1 = ii1(1); ii2 = ii2(end);
	if any(isnan(da(ii1,2:3))) && ~any(isnan(da(ii2,2:3)))
		dof = da(ii2,2:3);
		irf_log('proc','Delta off from the previous interval')
	elseif ~any(isnan(da(ii1,2:3))) && any(isnan(da(ii2,2:3)))
		dof = da(ii1,2:3);
		irf_log('proc','Delta off from the next interval')
	elseif any(isnan(da(ii1,2:3))) && any(isnan(da(ii2,2:3)))
		dof = [];
		irf_log('proc','No Delta offsets at this time')
	else
		da = da(~isnan(da(:,2)),:);
		dof = interp1(da(:,1),da(:,2:3),t,'spline');
	end
else
	dof = interp1(da(:,1),da(:,2:3),data(:,1),'spline');
	dof = [data(:,1) dof];
end


