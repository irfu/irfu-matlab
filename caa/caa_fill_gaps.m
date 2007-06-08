function res = caa_fill_gaps(data,te)
%CAA_FILL_GAPS(data,te)  fill gaps in the of a dataset
%
% res = caa_fill_gaps(data,te)
%       append NaNs to the end of DATA untill TE
%
% $Id$

% Copyright 2006-2007 Yuri Khotyaintsev

res = data;

if size(data,1)<2
	irf_log('proc','cannot fill gaps (not enough samples)')
	return
end

fs = c_efw_fsample(data);
if fs<=0
	irf_log('proc','cannot fill gaps (no sampling frequency)')
	return
end

ngap = fix((te -data(end,1))*fs);
if ngap>0
	tt = zeros(ngap+1,size(data,2));
	tt(:,1) = linspace(data(end,1),data(end,1)+ngap/fs,ngap+1);
	tt = tt(2:end,:);
	tt(:,2:end) = NaN;
	res = [res; tt];
	% Do not report 1 point gaps
	if ngap>1
		irf_log('proc',...
			sprintf('filling %d gaps at %s',ngap,epoch2iso(data(end,1),1)))
	end
end
