function [ii1,ii2]=irf_find_comm_idx(d1,d2)
%IRF_FIND_COMM_IDX   Find common index of two time series
%
% [ii1,ii2]=irf_find_comm_idx(e1,e2)
% 
% This function can be used when combining two datasets with
% gaps, like p12 and p34.
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev

DEBUG = 0;

ii1 = [];
ii2 = [];

% we use only time
t1 = d1(:,1);
t2 = d2(:,1);
s1 = 1; s2 = 1; % start index

n_loop = 0;

while 1
	la = min(length(t1(s1:end)),length(t2(s2:end))) - 1;
	d = t1(s1:s1+la)-t2(s2:s2+la);
	ii = find(abs(d) > 0);

	if isempty(ii) % time lines are identical
		% save the interval and break
		ii1 = [ii1 s1:s1+la];
		ii2 = [ii2 s2:s2+la];
		if DEBUG, disp('break at 1'), end
		break 
	end

	ii = ii(1); % first point where times differ
	if ii>1
		% save the interval and shift
		if DEBUG, disp('save interval and shift'), end
		ii1 = [ii1 s1:s1+ii-2];
		ii2 = [ii2 s2:s2+ii-2];
		s1 = s1 + ii - 1;
		s2 = s2 + ii - 1;
	end
	if d(ii) > 0
		% shift t2
		if DEBUG, disp('shift t2'), end

		i_t0 = find(t2(s2:end)==t1(s1));
		if isempty(i_t0)
			if t1(s1) > t2(end)
				if DEBUG, disp('break at 2'), end
				break
			end
			s1 = s1 + 1;
			it = find(t2>=t1(s1));
			s2 = it(1);
			clear it
		else
			s2 = s2 + i_t0 - 1;
		end
	else
		%shift t1
		if DEBUG, disp('shift t1'), end
		i_t0 = find(t1(s1:end)==t2(s2));
		if isempty(i_t0)
			if t2(s2) > t1(end)
				if DEBUG, disp('break at 3'), end
				break
			end
			s2 = s2 + 1;
			it = find(t1>=t2(s2));
			s1 = it(1);
			clear it
		else
			s1 = s1 + i_t0 - 1;
		end
	end
	n_loop = n_loop + 1;
	irf_log('proc',sprintf('gap # %d',n_loop))
end
