function phase_out = c_phase(t,phase_2)
%C_PHASE  Find spacecraft phase for given time vector
%
% phase_out = c_phase(t,phase_2)
%
% Find spacecraft phase for give time vector t
%
% t - column vector with time in isdat epoch
% phase_2 - column vector [time phase_2]
% phase_out - column vector [t phase]
%
% $Id$

error(nargchk(2,2,nargin))

if size(t,1)>1 && size(t,2)>1, error('t must be a vector'), end
if size(phase_2,1)<2, error('not enough points in phase_2'), end

MAX_SPIN_PERIOD = 4.3; % Sane value for Cluster

t=t(:); % t should be column vector
phase_out = [];

% interp1q sets points which are outside the time interval to NaN.
% We must check for this.
tt = t(t>=phase_2(1,1) & t<=phase_2(end,1));
if isempty(tt), return, end

while size(phase_2,1)>2
	ii = find(diff(phase_2(:,1))>MAX_SPIN_PERIOD);
	if isempty(ii)
		kk = find(tt>=phase_2(1,1) & tt<=phase_2(end,1));
		if ~isempty(kk)
			phase_out = [phase_out; tt(kk) interp1q(phase_2(:,1),phase_2(:,2),tt(kk))];
		end
		break
	else
		irf_log('proc',['gap in phase at ' epoch2iso(phase_2(ii(1),1),1)])
		if ii(1)==1
			irf_log('proc','throwing away 1 bad point')
			phase_2(1,:) = [];
		else
			kk = find(tt>=phase_2(1,1) & tt<=phase_2(ii(1),1));
			if ~isempty(kk)
				phase_out = [phase_out; ...
					tt(kk) interp1q(phase_2(1:ii(1),1),phase_2(1:ii(1),2),tt(kk))];
			end
			phase_2(1:ii(1),:) = [];
		end
	end
end

n_out = length(t) - length(phase_out);
if n_out > 0
	irf_log('proc',['throwing away ' num2str(n_out) ' points'])
end


