function out = SPDFTT2000toDATENUM(tt2000)
%SPDFTT2000toDATENUM converts the time in UTC string (returned from spdfcdfread)
%                or date/time values for CDF_TIME_TT2000 to MATLAB datenum
%
%   OUT = SPDFTT2000toDATENUM(tt2000) returns MATLAB datenum.
%   OUT a column vector of numerical values of MATLAB date numbers.
%
%     tt2000               A vector or cell of UTC string or an either
%                          M-by-1 matrix of CDF_TIME_TT2000 values or
%                          M-by-9 matrix containing M rows, each with date/time
%                          fields for year, month, day, hour, minute, second,
%                          millisecond, microsecond, and nanosecond, in that
%                          order.
%
%   Note:
%     The valid tt2000 string should have one of the following forms:
%     1. dd-mmm-yyyy hh:mm:ss.mmmuuunnn (length of 30), e.g.,
%       "01-JAN-2000 12:00:00.123456789"
%     2. yyyymmdd.dddddddddd (length of 19), e.g.,
%       "20000101.1200000000"
%     3. yyyymmddhhmmss (length of 14), e.g.,
%       "20000101120000"
%     4. yyyy-mm-ddThh:mm:ss.mmmuuunnn (ISO 8601, length of 29), e.g.,
%       "2000-01-01T12:00:00.123456789"
%     where mmmuuunnn is milliseconds, microseconds and nanoseconds.
%
%   Examples:
%
%   % Read all the variable data in a CDF file. Among them, the variable of 
%     CDF_TIME_TT2000 data type, at index of 17, is returned. Convert the
%     variable of CDF_TIME_TT2000 data type in cdftt2000 objects to MATLAB
%     datenum.
%
%   data = spdfcdfread('test','KeepEpochAsIs',true);
%   tt2000 = data(1,17);
%   datenums = SPDFTT2000toDATENUM(tt2000);
%
%   % Convert the UTC strings in vector to MATLAB datenum.
%
%   tt2000 = ['2009-01-01T00:00:00.123456789';
%             '2009-01-01T12:00:00.123456789'];
%   datenums = SPDFTT2000toDATENUM(tt2000);
%
%   % Convert the date/times in matrix to MATLAB datenum.
%
%   tt2000 = [2009 01 01 00 00 00 123 456 789;
%             2009 01 01 12 00 00 123 456 789];
%   datenums = SPDFTT2000toDATENUM(tt2000);
%
%   See also CDFTT2000, SPDFENCODETT2000, SPDFCOMPUTETT2000, SPDFPARSETT2000

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFTT2000toDATENUM:inputArgumentCount', ...
          'SPDFTT2000toDATENUM requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFTT2000toDATENUM:outputArguments', ...
          'SPDFTT2000toDATENUM requires only one output argument.')
end

if (isa(tt2000, 'cdftt2000'))
  error('MATLAB:SPDFTT2000toDATENUM:inputobject', ...
        'SPDFTT2000toDATENUM can not handle input data in cdftt2000 object class.')
end

s = size(tt2000);
if (iscell(tt2000(1,1)))
  if (s(1) == 1)
    tt2000 = tt2000{1,1};
  else
    if (iscellstr(tt2000))
      tt2000 = char(tt2000);
    else
      if (s(2) == 1)
        for x=1:s(1)
        tt2000(x,1) = tt2000{x,1};
        end
      end
    end
  end
end
if (iscell(tt2000))
  tt2000 = cell2mat(tt2000);
end
if (isnumeric(tt2000))
  s = size(tt2000);
  if (s(2) == 1)
    dates = spdfbreakdowntt2000(tt2000);
    dates3 = dates(:, 1:6);
    for p =1:s(1)
      dates3(p,6) = dates(p,6) + dates(p,7)/1000 + dates(p,8)/1000000 + ...
                    dates(p,9)/1000000000;
    end
    out = datenum(dates3);
  else
    if (s(2) ~= 9)
      error('MATLAB:SPDFTT2000toDATENUM:inputArgumentvector', ...
            'SPDFTT2000toDATENUM requires an M by 9 vector for date/time fields.')
    end
    dates=tt2000(:,1:6);
    for p =1:s(1)
      dates(p,6) = tt2000(p,6) + tt2000(p,7)/1000 + tt2000(p,8)/1000000 + ...
                   tt2000(p,9)/1000000000;
    end
    out = datenum(dates); 
  end
else
  dates = spdfparsett2000(tt2000);
  dates2 = spdfbreakdowntt2000(dates);
  dates3 = dates2(:, 1:6);
  for p =1:s(1)
    dates3(p,6) = dates2(p,6) + dates2(p,7)/1000 + dates2(p,8)/1000000 + ...
                  dates2(p,9)/1000000000;
  end
  out = datenum(dates3);
end
