%
% Convert from MJS time to date 
%
%
function [year,mon,day,hour,min,sec]=cefMjsToDate(mjs)

  temp = mjs/86400 + 0.5/86400000
  jday = floor(temp);
  l = floor(floor(4000*(jday  + 18204))/1461001);
  n = jday - floor(floor(1461*l)/4) + 18234;
  m = floor(floor(80*n)/2447);
  day = n - floor(floor(2447*m)/80);
  jj = floor(m/11);
  mon = m + 2 -12*jj;
  year = 1900 + l + jj;
  temp = (temp-jday)*24
  hour =  floor(temp)
  temp = (temp -(hour))*60;
  min =  floor(temp)
  temp = (temp -(min))*60.0;
  sec =  temp

end