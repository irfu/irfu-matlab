function [ii1,ii2]=findCommInd(d1,d2)
%findCommInd find common index of two time series
% [ii1,ii2]=findCommInd(e1,e2)
% 
% This function can be used when combining two datasets with
% gaps, like p12 and p34.
%
% $Revision$  $Date$
%

% Copyright 2004 Yuri Khotyaintsev

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

	if isempty(ii)
		% save the interval and break
		ii1 = [ii1 s1:s1+la];
		ii2 = [ii2 s2:s2+la];
		break 
	end

	ii = ii(1); % first point where times differ
	if ii>1
		% save the interval and shift
		ii1 = [ii1 s1:s1+ii-2];
		ii2 = [ii2 s2:s2+ii-2];
		s1 = s1 + ii - 1;
		s2 = s2 + ii - 1;
	end
	if d(ii) > 0
		% shift t2
		i_t0 = find(t2(s2:end)==t1(s1));
		if isempty(i_t0), break, end
		s2 = s2 + i_t0 - 1;
	else
		%shift t1
		i_t0 = find(t1(s1:end)==t2(s2));
		if isempty(i_t0), break, end
		s1 = s1 + i_t0 - 1;
	end
	n_loop = n_loop + 1;
	disp(sprintf('gap # %d',n_loop))
end
