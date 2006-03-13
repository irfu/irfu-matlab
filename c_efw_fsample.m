function sf = c_efw_fsample(data)
%C_EFW_FSAMPLE  guess sampling frequency
%
% fs = c_efw_fsample(data)
% returns sampling frequency in Hz
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

data = data(:,1);

if length(data)<=1
	error('cannot compute sampling frequency from less than two points')
end

sf = guess_fsample((length(data) - 1)/(data(end) - data(1)));
if sf, return
else
	sf = guess_fsample(1/(data(2) - data(1)));
	if sf, return
	elseif length(data)>2
		sf = guess_fsample(1/(data(3) - data(2)));
		if sf, return
		else
			sf = guess_fsample(1/(data(end) - data(end-1)));
			if sf, return, end
		end
	end
end

if ~sf, irf_log('proc','cannot guess sampling frequency'), end
return

function sf = guess_fsample(f)
K_PLUS = 1.1;
K_MINUS = .9;

if f<K_PLUS*5 & f>K_MINUS*5, sf = 5;          % LX
elseif f<K_PLUS*25 & f>K_MINUS*25, sf = 25;      % NM
elseif f<K_PLUS*450 & f>K_MINUS*450, sf = 450;   % BM1
elseif f<K_PLUS*9000 & f>K_MINUS*9000, sf = 9000;% IB
elseif f<K_PLUS*.25 & f>K_MINUS*.25, sf = .25;   % SPIN
else, sf = 0;
end
