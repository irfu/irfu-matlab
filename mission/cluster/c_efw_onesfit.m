function spinfit = c_efw_onesfit(pair,fout,maxit,minpts,te,data,tp,ph)
%C_EFW_ONESFIT produce one spin fit value (ex, ey)
%   of EFW data frome given probe pair. The spin fit will
%   be to all the data put in, which thus typically should
%   cover exactly one spin.
%
% spinfit = c_efw_onesfit(pair,fout,maxit,minpts,te,data,tp,ph)
%
% Input:
%  pair - probe pair used (12, 32, 34)
%  fout - minimum fraction of fit standard deviation that defines an outlier
%         (zero means no removal of outliers)
%  maxit - maximum number of iterations (zero means infinity)
%  minpts - minimum number of data points to perform fit
%      (set to 5 if smaller number if given)
%  te - EFW time in seconds (isGetDataLite time)
%  data - EFW data from pair in mV/m, should correspond to te
%  tp - Ephemeris time in seconds (isGetDataLite time)
%  ph - Ephemeris phase in degr (sun angle for s/c Y axis), should
%      correspond to tp
%
% Output:
%  spinfit = [ts,ex,ey,offset,sdev0,sdev,iter,nout]
%  ts - time vector in seconds
%  ex - E-field x-component in DSI coordinates (almost GSE)
%  ey - E-field y-component in DSI coordinates (almost GSE)
%  offset - mean value of input data
%  sdev0 - standard deviation in first fit
%  sdev - standard deviation in final fit
%  iter - number of iterations (one if OK at once)
%  nout - number of outliers removed
% If fit was unsuccesful, all output except iter is zero.
%
% See also C_EFW_SFIT
%

% Anders.Eriksson@irfu.se, 13 December 2002

% Defaults:
if(minpts < 5)
  minpts = 5;
end

% Make columns of all input data:
te = tocolumn(te);
data = tocolumn(data);
tp = tocolumn(tp);
ph = tocolumn(ph);

% Calcluate phase (in rad) at EFW sample times:
ph = unwrap(pi*ph/180);
pol = polyfit(tp,ph,1);
ph = polyval(pol,te);

% Find phase of given pair:
pair_ok = 1;
if pair == 12
  ph = ph + 3*pi/4;
elseif pair == 32
  ph = ph + pi/2;
elseif pair == 34
  ph = ph + pi/4;
else
  pair_ok = 0;
end

% Do spin fit

spinfit = zeros(1,8);
ph1 = ph;
data1 = data;
if pair_ok
  ready = 0;
  iter = 0;
  nout = 0;
  while(~ready)
    iter = iter + 1;
    if max(size(data1)) < minpts
      spinfit(7) = iter;
      ready = 1;
    else
      spfit = regress(data1, [cos(ph1) sin(ph1) ones(size(ph1))]);
      fit = spfit(1) * cos(ph1) + spfit(2) * sin(ph1) + spfit(3);
      dif = data1 - fit;
      sdev = std(dif);
      if iter == 1, sdev0 = sdev; end
      ind = find(abs(dif) > fout * sdev);  % Find outliers
      if fout == 0 || maxit == 0 || isempty(ind) || iter >= maxit
        spinfit(1) = mean(te);
        spinfit(2) = spfit(1);
        spinfit(3) = -spfit(2);  % Because s/c is spinning upside down
        spinfit(4) = spfit(3);
        spinfit(5) = sdev0;
        spinfit(6) = sdev;
        spinfit(7) = iter;
        spinfit(8) = nout;
        ready = 1;
      else
        ph1(ind) = [];
        data1(ind) = [];
        nout = nout + max(size(ind));
      end
    end
  end
end
