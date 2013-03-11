function out = UTC2CDFEpoch(UTC,month,day,hour,minute,second,millisec)
%UTC2CDFEpoch converts a UTC date/time in string or components to CDF_EPOCH
%
%   There are two forms of this function:
%
%   OUT = UTC2CDFEpoch(UTC) 
%         Parses a single CDF epoch string and returns a CDF_EPOCH
%         data type.
%
%     UTC               A UTC string
%
%   Note: the valid epoch string should be one of the following forms:
%      0: dd-mmm-yyyy hh:mm:ss.lll, e.g., "01-JAN-2010 12:00:00.000"
%      1: yyyymmdd.ddddddd, e.g., "20100101.1200000"
%      2: yyyymmddhhmmss, e.g., "20100101120000"
%      3: yyyy-mm-ddThh:mm:ss.lll, e.g., "2010-01-01T12:00:00.000"
%         where lll is milliseconds.
%   Or,
%
%   OUT = UTC2CDFEpoch(year,month,day,hour,minute,second,millisec) 
%         Compute the CDF epoch in CDF_EPOCH data type.
%
%     year, month, day, hour,       Integer form for date/time components
%     minute, second, millisecond,  
%     microsecond, nanosecond
%
%   Examples:
%
%   % Write 100 epoch records, starting from 2009-01-01T00:00:00.123 with
%   % one (1) second stride, to a CDF.
%
%   utc = '2009-01-01T00:00:00.123';
%   epoch = UTC2CDFEpoch(utc);
%   epochs = epoch+[0:99]*1000;
%   CDFWRITE('example', {'Epoch', epochs}, 'EPOCHArrayisCDFEpoch', true, ...
%            'recordbound', {'Epoch'});
%
%   % Alternatively, the sample can be entered as
%
%   epoch = UTC2CDFEpoch(2010,1,1,0,0,0,0);
%   epochs = epoch:[0:99]*1000;
%   CDFWRITE('example', {'Epoch', epochs}, 'EPOCHArrayisCDFEpoch', true, ...
%            'recordbound', {'Epoch'});
%
%   See also CDFWRITE, CDFREAD, CDFEPOCH.

%   Copyright 1984-2009 The MathWorks, Inc.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

function day = JulianDay (y,m,d)
  t1=floor((m+9)/12);
  t2=floor((7*(y+t1))/4);
  if ((m-9) > 0)
    l1=floor((m-9)/7);
  else
    l1=ceil((m-9)/7);
  end
  l2=floor((y+l1)/100)+1;
  l3=floor((3*l2)/4);
  y1=floor((275*m)/9);
  day=367*y-t2-l3+y1+d+1721029;
end

function mon = MonthToken (month)
  tmp= lower(month);
  if (isequal(tmp,'jan')) mon = 1;
  elseif (isequal(tmp,'feb')) mon = 2;
  elseif (isequal(tmp,'mar')) mon = 3;
  elseif (isequal(tmp,'apr')) mon = 4;
  elseif (isequal(tmp,'may')) mon = 5;
  elseif (isequal(tmp,'jun')) mon = 6;
  elseif (isequal(tmp,'jul')) mon = 7;
  elseif (isequal(tmp,'aug')) mon = 8;
  elseif (isequal(tmp,'sep')) mon = 9;
  elseif (isequal(tmp,'oct')) mon = 10;
  elseif (isequal(tmp,'nov')) mon = 11;
  elseif (isequal(tmp,'dec')) mon = 12;
  else mon = 0;
  end
end

function epoch = computeepoch (yy,mm,dd,hh,nn,ss,ll)
  if (yy==9999 && mm==12 && dd==31 && hh==23 && nn==59 && ss==59 && ...
      ll==999)
    epoch = -1.0E31;
  else
    daysfrom0AD = JulianDay(yy,mm,dd) - 1721060;
    epoch = daysfrom0AD * 86400000.0 + hh * 3600000.0 + nn * 60000.0 + ...
          ss * 1000.0 + ll;
  end
end

if (nargin ~= 1 && nargin ~= 7)
    error('MATLAB:UTC2CDFEpoch:inputArgumentCount', ...
          'UTC2CDFEpoch requires one input string or seven date/time components.')
end

if (nargout > 1)
    error('MATLAB:UTC2CDFEpoch:outputArguments', ...
          'UTC2CDFEpoch requires only one output argument.')
end

if (nargin == 1 && ischar(UTC))
  if (length(UTC) == 24) 
    yy = str2num(UTC(8:11));
    mon = UTC(4:6);
    mm=MonthToken(mon);
    dd = str2num(UTC(1:2));
    hh = str2num(UTC(13:14));
    nn = str2num(UTC(16:17));
    ss = str2num(UTC(19:20));
    ll = str2num(UTC(22:24));
    out = computeepoch(yy,mm,dd,hh,nn,ss,ll);
  elseif (length(UTC) == 16) 
    yy = str2num(UTC(1:4));
    mm = str2num(UTC(5:6));
    dd = str2num(UTC(7:8));
    hh = str2num(UTC(10:end));
    len =length(UTC(10:end));
    ll = 86400000*(hh/(10^(len+1)));
    out = computeepoch(yy,mm,dd,0,0,0,ll);
  elseif (length(UTC)  == 14)
    yy = str2num(UTC(1:4));
    mm = str2num(UTC(5:6));
    dd = str2num(UTC(7:8));
    hh = str2num(UTC(9:10));
    nn = str2num(UTC(11:12));
    ss = str2num(UTC(13:14));
    out = computeepoch(yy,mm,dd,hh,nn,ss,0);
  elseif (length(UTC) == 23)
    yy = str2num(UTC(1:4));
    mm = str2num(UTC(6:7));
    dd = str2num(UTC(9:10));
    hh = str2num(UTC(12:13));
    nn = str2num(UTC(15:16));
    ss = str2num(UTC(18:19));
    ll = str2num(UTC(21:23));
    out = computeepoch(yy,mm,dd,hh,nn,ss,ll);
  else
    error('MATLAB:UTC2CDFEpoch:string', ...
          'UTC2CDFEpoch string does not meet requirements.')
  end
elseif (nargin == 7) 
  yy = UTC;
  mm = month;
  dd = day;
  hh = hour;
  nn = minute;
  ss = second;
  ll = millisec;
  out = computeepoch(yy,mm,dd,hh,nn,ss,ll);
else
  error('MATLAB:UTC2CDFEpoch:input', ...
        'UTC2CDFEpoch input does not meet requirements.')
end
end
