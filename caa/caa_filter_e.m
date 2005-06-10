function res = caa_filter_e(data,wind)
% CAA_FILTER_E  simple high pass filter EFW E
%
% res = caa_filter_e(data,[window])
%	simple high pass filter, will removemoving average over WINDOW
%	WINDOW in seconds, delault 0.5 sec
%
% $Id$

if nargin<2, wind=.5; end

ii = find( ~isnan(data(:,2)) );

sf = c_efw_fsample(data(ii,1));
nw2 = ceil(sf*wind/2);

% Order data
ndata = round((data(end,1)-data(1,1))*sf+1);
t2 = data(1,1) + (ndata-1)/sf;

nkomp = size(data,2)-1;
value = 1e10;
E = ones(ndata,nkomp)*value;

ind = round((data(:,1)-data(1,1))*sf+1);
E(ind,:) = data(:,2:end);

E(find(E(:,1)==value),:) = NaN;
ttt = ones(nw2,nkomp)*NaN;
E_tmp = [ttt; E; ttt];

EE = zeros(ndata,2*nw2+1);

for komp = 1:nkomp
	for k=1:ndata, EE(k,:) = E_tmp(k:k+nw2*2,komp); end
	m = mean(EE,2);
	E(:,komp) = E(:,komp) - m;
end

res = data;
res(:,2:end) = E(ind,:);
