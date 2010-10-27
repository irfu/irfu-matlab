function [ t0, dt, tplus, tminus, corr ] = mycorr( data1, data2, fs )
%MYCORR  Cross correlation "by eye"
%
%   Minimize S = SUM ( [data1 - data2]^2 )
%
%   [ t0, dt, tplus, tminus, corr ] = mycorr( data1, data2, fs )
%
%   If x and y are not the same length, the shorter vector is zero-padded
%   to the length of the longer vector
%
%   t0 - time delay between data1 and data2
%   dt = (tplus + tminus)/2

if nargin <3, fs = 1; end

ndata = max(size(data1,1),size(data2,1));

if ndata < 5, error('Must have at least 5 points!'), end
    
ncomp = size(data1,2);

% Remove mean - XXX: seems to increase the error
%data1 = data1 - ones(size(data1,1),1)*mean(data1,1);
%data2 = data2 - ones(size(data2,1),1)*mean(data2,1);

if length(data1) < ndata
    data1(end:ndata,:) = NaN;
elseif length(data2) < ndata
    data2(end:ndata,:) = NaN;
end

corr = zeros(ndata*2-1,1);

data2 = [NaN*ones(ndata-1,ncomp); data2; NaN*ones(ndata-1,ncomp)];


for i=(-ndata+1):ndata-1
    S = ( data1 - data2(i+ndata:i+2*ndata-1,:) ).^2;
    if any(~isnan(S))
        corr(i+ndata) = sum(S(~isnan(S)))/(ndata - abs(i))/ncomp;
    else corr(i+ndata) = NaN;
    end
end

OFF = 1;
ii = find(corr == min(corr(1+OFF:end-OFF))); % take out the last points for security

warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
p = polyfit((ii-OFF:ii+OFF)',corr(ii-OFF:ii+OFF),2);
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
t0 = -p(2)/p(1)/2;
if abs(t0-ii)<1
    corr_min = polyval(p,t0);
else
    disp('Discarding parabolic fit!')
    t0 = ii;
    corr_min = corr(ii(1));
end

% We define the error as DT at which the correlation doubles
corr_plus = corr - 2*corr_min;

ii_plus = find(corr_plus > 0);
if isempty(ii_plus), tplus = NaN; tminus = NaN;
else
    ii_plus = ii_plus - ii;
    ii1 = find(ii_plus>0);
    ii1 = ii_plus(ii1(1)) + ii;
    tplus = find_zero(corr_plus,ii1,ii1-1);
    ii1 = find(ii_plus<0);
    ii1 = ii_plus(ii1(end)) + ii;
    tminus = find_zero(corr_plus,ii1,ii1+1);
end

dt = ( tplus - tminus) /2;
tplus = tplus - ndata;
tminus = tminus - ndata;
t0 = t0 - ndata;

% Diagnostics
if false
    figure; %#ok<UNRCH>
    dt_tmp = (1:(ndata*2-1)) - ndata;
    plot(dt_tmp/fs,corr,'.-')
    set(gca, 'YLim',[0 3*corr_min])
    hold on
    plot(t0/fs,corr_min,'r*')
    text(t0/fs,corr_min,sprintf(' t0 = %.3f s +/- %.3f s (%.2f +/- %.2f)',t0/fs,dt/fs,t0,dt))
    plot(tplus/fs,corr_min*2,'b*')
    text(tplus/fs,corr_min*2,sprintf(' t+ = %.2f s (%.2f)',tplus/fs,tplus))
    plot(tminus/fs,corr_min*2,'b*')
    text(tminus/fs,corr_min*2,sprintf(' t- = %.2f s (%.2f)',tminus/fs,tminus))
    hold off
    grid
end

% Prepare the output
tplus = - tplus/fs;
tminus = - tminus/fs;
t0 = - t0/fs;
dt = dt/fs;

end
