function epochSec = toepoch(x)
% toepoch - Convert a [YYYY MM DD hh mm ss] time specification to seconds
% since 1970/01/01.
%
% If YYYY is less than 100 it is assumed to bethe 20th century.
% (ie [80 1 1 0 0 0] => [1980 1 1 0 0 0], in other words: Jan 1st, 1980.
%
% If input is less than the full vector it will be padded with zero, for
% seconds, minutes and hours. (ie. [1980 1 1] => [1980 1 1 0 0 0], in other
% words: midnight Jan 1st, 1980.
%
% Note: Leap seconds are NOT taken into account. 
% 

[m, n] = size(x);
if n~=2 && n~=3 && n~=6,
  if m==2 || m==3 || n==6
    irf.log('warning','Unexpected vector input in toepoch. Be careful with results.');
     x = x';
     n = size(x,2);
  else
    irf.log('warning','Illegal argument.');
    epochSec = NaN;
    return
  end
end

if n==6
  if x(1,1)<100
    irf.log('warning','Assuming 20th century. Be careful with results.');
    x(:,1) = 1900 + x(:,1);
  end
end

% Unix epoch is actually to be defined as a signed "int32" but alot of our
% code assumes only double for unix epoch.
epochSec = round(86400 *( datenum(x) - datenum('01-Jan-1970')));

end