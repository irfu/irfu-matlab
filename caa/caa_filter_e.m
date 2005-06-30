function res = caa_filter_e(data,wind)
% CAA_FILTER_E  simple high pass filter EFW E
%
% res = caa_filter_e(data,[window])
%	simple high pass filter, will removemoving average over WINDOW
%	WINDOW in seconds, delault 0.5 sec
%
% $Id$

MAXDATA = 100000;

if nargin<2, wind=.5; end

ii = find( ~isnan(data(:,2)) );

sf = c_efw_fsample(data(ii,1));
nw2 = ceil(sf*wind/2);

% Order data
ndata = round((data(end,1)-data(1,1))*sf+1);

nkomp = size(data,2)-1;
value = 1e10;
E = ones(ndata,nkomp)*value;

ind = round((data(:,1)-data(1,1))*sf+1);
E(ind,:) = data(:,2:end);

E(find(E(:,1)==value),:) = NaN;
ttt = ones(nw2,nkomp)*NaN;
E_tmp = [ttt; E; ttt];
clear ttt

n_start = 1;
if ndata>MAXDATA; n_end = MAXDATA;
else, n_end = ndata;
end

EE = zeros(n_end,2*nw2+1);
while n_start<=ndata
	for komp = 1:nkomp
		for k=n_start:n_end, EE(k-n_start+1,:) = E_tmp(k:k+nw2*2,komp); end
		m = mean(EE(1:n_end-n_start+1,:),2);
		E(n_start:n_end,komp) = E(n_start:n_end,komp) - m;
	end
	if n_end==ndata, break, end
	n_start = n_end + 1;
	n_end = n_end + MAXDATA;
	if n_end>ndata, n_end = ndata; end
end
clear EE

res = data;
res(:,2:end) = E(ind,:);
