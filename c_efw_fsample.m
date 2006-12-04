function sf = c_efw_fsample(data,mode)
%C_EFW_FSAMPLE  guess sampling frequency
%
% fs = c_efw_fsample(data,[mode])
% returns sampling frequency in Hz
% mode - restrict to particular type of data :
%        'any' - default 0.25/5/25/450/4500/9000 Hz
%        'lx'  - LX data 5 Hz
%        'hx'  - HX data 25/450 Hz
%        'ib'  - internal burst data 450/4500/9000 Hz
%
% $Id$

% Copyright 2005,2006 Yuri Khotyaintsev

data = data(:,1);

if length(data)<=1
	error('cannot compute sampling frequency from less than two points')
end

if nargin <2, mode='any'; end
if ~(strcmp(mode,'any') || strcmp(mode,'lx') || strcmp(mode,'hx') || strcmp(mode,'ib'))
	error('bad value for MODE')
end

sf = guess_fsample((length(data) - 1)/(data(end) - data(1)),mode);
if sf, return
else
	sf = guess_fsample(1/(data(2) - data(1)),mode);
	if sf, return
	elseif length(data)>2
		sf = guess_fsample(1/(data(3) - data(2)),mode);
		if sf, return
		else
			sf = guess_fsample(1/(data(end) - data(end-1)),mode);
			if sf, return, end
		end
	end
end

if ~sf, irf_log('proc','cannot guess sampling frequency'), end
return

function sf = guess_fsample(f,mode)
K_PLUS = 1.1;
K_MINUS = .9;

if f<K_PLUS*5 && f>K_MINUS*5 && (strcmp(mode,'any') || strcmp(mode,'lx'))
	sf = 5;     % LX
elseif f<K_PLUS*25 && f>K_MINUS*25 && (strcmp(mode,'any') || strcmp(mode,'hx'))
	sf = 25;    % NM
elseif f<K_PLUS*450 && f>K_MINUS*450 && (strcmp(mode,'any') || strcmp(mode,'hx') || strcmp(mode,'ib'))
	sf = 450;   % BM1/IB
elseif f<K_PLUS*4500 && f>K_MINUS*4500 && (strcmp(mode,'any') || strcmp(mode,'ib'))
	sf = 4500;  % IB
elseif f<K_PLUS*9000 && f>K_MINUS*9000 && (strcmp(mode,'any') || strcmp(mode,'ib'))
	sf = 9000;  % IB
elseif f<K_PLUS*.25 && f>K_MINUS*.25 && strcmp(mode,'any'), sf = .25;   % SPIN
else sf = 0;
end
