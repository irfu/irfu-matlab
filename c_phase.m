function phase_out = c_phase(t,phase_2)
%C_PHASE spacecraft phase for give time vector
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

if size(t,1)>1 & size(t,2)>1, error('t must be a vector'), end

t=t(:); % t should be column vector

% interp1q sets points which are outside the time interval to NaN.
% We must check for this.
ii = find(t>=phase_2(1,1) & t<=phase_2(end,1));
if ~isempty(ii)
	phase_out = [t(ii) interp1q(phase_2(:,1),phase_2(:,2),t(ii))];
	n_out = length(t) - ii(end) + ii(1) - 1;
	if n_out > 0
		irf_log('proc',['throwing away ' num2str(n_out) ' points'])
	end
else 
	phase_out = [];
end

