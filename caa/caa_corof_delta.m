function data = caa_corof_delta(data,probe_p,dof,action)
%CAA_COROF_DELTA correct delta offsets
%
% new_data = caa_corof_delta(data,deltaoff)
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if probe_p~=12 && probe_p~=32 && probe_p~=34
	error('bag value for PROBE_P')
end

if probe_p==34 && isreal(dof), return, end
if (probe_p == 12 || probe_p == 32) && ~isreal(dof), return, end

if nargin<4, action = 1;
else
	if ischar(action)
		if strcmp(action,'apply'), action = 1;
		elseif strcmp(action,'undo'), action = -1;
		end
	elseif action~=1 || action~=-1
		error('bag value for ACTION')
	end
end

% Real offset is applied to p12/32, imaginary to p34.
for comp=1:2
	if (isreal(dof(comp)) && probe_p~=34) || (~isreal(dof(comp)) && probe_p==34)
		if ~isreal(dof(comp)), dof(comp) = imag(dof(comp)); end
		data(:,comp) = data(:,comp) - action*dof(comp);
	end
end


