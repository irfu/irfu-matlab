function res=irf_find_max(data, fla)
%IRF_FIND_MAX  Find maxima in data
%
% res=irf_find_max(data, [flag])
% Find maxima using quadratic and cubic interpolation.
% If flag=1, minima are found instead.
%

% Part of the code was copied from /usr/local/matl7/toolbox/optim/findmax.m

[n,m]=size(data);
if m>2, error('input must be 1D: [time data]'),end

if nargin<2, fla=0; end

t = data(:,1);
if fla, data = -data(:,2);
else, data = data(:,2);
end

ii = find( (data>[data(1)-1;data((1:n-1)')]) & (data>=[data((2:n)');data(n)-1]) );
s = length(ii);
res_m = zeros(s,1);
res_t = zeros(s,1);

for i=1:s
  ix=ii(i);
  if (ix==1 || ix==n)
    res_m(i) = data(ix);
    res_t(i) = t(ix);
  elseif ix>n-3 || ix<4
    if data(ix-1)==data(ix) && data(ix+1)==data(ix)
      res_m(i) = data(ix);
      res_t(i) = t(ix);
    else
      [res_t(i), res_m(i)] = quad_max(t(ix-1:ix+1),data(ix-1:ix+1));
    end
  elseif data(ix+2)<data(ix+1)
    [res_t(i), res_m(i)] = cub_max(t(ix-1:ix+2),data(ix-1:ix+2));
  elseif data(ix-2)<data(ix-1)
    [res_t(i), res_m(i)] = cub_max(t(ix-2:ix+1),data(ix-2:ix+1));
  elseif data(ix)>data(ii(i+1))
    [res_t(i), res_m(i)] = quad_max(t(ix-1:ix+1),data(ix-1:ix+1));
  else
    res_m(i) = data(ix);
    res_t(i) = t(ix);
  end
end

if fla, res = [res_t -res_m];
else, res = [res_t res_m];
end

function [tm, m] = quad_max(t,d)

% Normalize t as described in help polyfit
t0 = mean(t);
tint = std(t);
t = (t - t0)/tint;

p = polyfit(t,d,2);
tm = -p(2)/p(1)/2;
m = p(1)*tm^2 + p(2)*tm + p(3);
tm = tm*tint + t0;

if tm<t(1) || tm>t(end), error('quad_max: max is outside the region'), end
%disp(['quad_max: ' num2str(tm) ' ' num2str(m)])

function [tm, m] = cub_max(t,d)

% Normalize t as described in help polyfit
t0 = mean(t);
tint = std(t);
t = (t - t0)/tint;

p = polyfit(t,d,3);
d = 4*p(2)^2 - 4*3*p(1)*p(3);
if d<0, error('no roots'), end
t1 = (-2*p(2) + sqrt(d))/p(1)/6;
t2 = (-2*p(2) - sqrt(d))/p(1)/6;

if (t1>=t(1) && t1<=t(end)), tm = t1;
elseif (t2>=t(1) && t2<=t(end)), tm = t2;
else
  error('quad_max: max is outside the region'),
end

m = p(1)*tm^3 + p(2)*tm^2 + p(3)*tm + p(4);
tm = tm*tint + t0;
%disp(['cub_max: ' num2str(tm) ' ' num2str(m)])
