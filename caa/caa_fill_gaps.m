function [data, filldata] = caa_fill_gaps(data,te)
%CAA_FILL_GAPS(data,te)  fill gaps in the of a dataset
%
% res = caa_fill_gaps(data,te)
%       append NaNs to the end of DATA untill TE
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,2,nargin))

if size(data,1)<2
	irf_log('proc','cannot fill gaps (not enough samples)')
	return
end

fs = c_efw_fsample(data);
if fs<=0
	irf_log('proc','cannot fill gaps')
	return
end

ngap = fix( (te -data(end,1))*fs - 1);
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
