function [secs]=toepoch(x)
% toepoch - Convert a [YYYY MM DD hh mm ss] time specification to seconds since 1970.
[m,n]=size(x);
if n~=2 && n~=3 && n~=6
  if m==2 || m==3 || n==6
     x=x'; n=size(x,2);
  else
    irf.log('warning','Illegal argument:\n')
   secs=NaN;
   return
  end
end

if n==2
  y=x;
  x(:,[3,6])=rem(y,100);
  x(:,[2,5])=rem(floor(y/100),100); 
  x(:,[1,4])=floor(y/10000);   
elseif n==3
  x(:,6)=rem(x(:,3),100); x(:,5)=floor(x(:,3)/100);
  x(:,4)=rem(x(:,2),100); x(:,3)=floor(x(:,2)/100);
  x(:,2)=rem(x(:,1),100); x(:,1)=floor(x(:,1)/100)+1900;
elseif n==6
  if x(1,1)<100, x(:,1)=1900+x(:,1); end
end 

years=x(:,1);

if isnan(years)
   secs=NaN;
   return
end

if(verLessThan('matlab','8.4'))
  for year=unique(x(:,1))'
    daym=[0 31 28 31 30 31 30 31 31 30 31 30 31];
    if leap_year(year), daym(3) = 29; end
    days = cumsum(daym(1:12))';
    ind = find(years==year);
    hours(ind,1) = (days(x(ind,2))+(x(ind,3)-1))*24+x(ind,4);
    secs(ind,1) = (hours(ind)*60+x(ind,5))*60+x(ind,6);
    yr_sec = 0;
    diff_yr = year - 1970;
    if(diff_yr>0)
      for i = 1:diff_yr
        if leap_year(1969+i)
            yr_sec = yr_sec + 31622400;
        else
            yr_sec = yr_sec + 31536000;
        end
      end
    elseif(diff_yr<0)
      for i = -1:-1:diff_yr
        if leap_year(1970+i)
          yr_sec = yr_sec - 31622400;
        else
          yr_sec = yr_sec - 31536000;
        end
      end
    end
    secs(ind,1) = secs(ind,1) + yr_sec;
  end
else
  % Use built in Matlab function (as of Matlab 8.4, R2014b)
  secs = posixtime(datetime(x));
end

% Help function
function leap = leap_year(year)
  leap = false(size(year)); % Assume no leap years
  indLeap = rem(year,4)==0;
  leap(indLeap) = true; % Assume these are leap years.
  indLeap = rem(year,100)==0 && rem(year,400)~=0;
  leap(indLeap) = false; % Not leap years: every even century not evenly divisable with 400, ie. 1700, 1800, 1900, 2100, 2200
end

end

