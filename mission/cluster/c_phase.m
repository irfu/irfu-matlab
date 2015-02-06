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

narginchk(2,2)

if size(t,1)>1 && size(t,2)>1, error('t must be a vector'), end
if size(phase_2,1)<2, error('not enough points in phase_2'), end

MAX_SPIN_PERIOD = 4.3; % Sane value for Cluster
MIN_SPIN_PERIOD = 3.6;

t=t(:); % t should be column vector
phase_out = [];

% Sanity check
badt=find(phase_2(:,1)<iso2epoch('2000-01-01T00:00:00Z'));
if ~isempty(badt)
	irf_log('proc',['Bad time ignored ' epoch2iso(phase_2(badt(1),1))])
	phase_2(badt,:) = [];
	clear badt;
end

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

% Sanity check
if ~isempty(phase_out())
    phase_unwrapped = unwrap(phase_out(:,2)/180*pi);
    SpinRate = diff(phase_unwrapped)./diff(phase_out(:,1));
    ii = find( diff(phase_out(:,1))< 0.95*2*pi/median(SpinRate) &...
        (SpinRate<2*pi/MAX_SPIN_PERIOD | SpinRate>2*pi/MIN_SPIN_PERIOD) );
end

if ~isempty(ii)
    ii_jump = find(diff(ii')>1);
    ii_jump = [1 ii_jump];
    for i=1:length(ii_jump)
        if i==length(ii_jump)
            kk = [ii(ii_jump(i)) ii(end)];
        else
            kk = [ii(ii_jump(i)) ii(ii_jump(i+1))-1];
        end
        if kk(1,2)==0, kk(1,2)=1; end;
        irf_log('proc',['bad phase at ' irf_disp_iso_range(phase_out(kk,1)',1)])
    end
    phase_out(ii,:) = [];
end

n_out = length(t) - length(phase_out);
if n_out > 0
	irf_log('proc',['throwing away ' num2str(n_out) ' points'])
end


