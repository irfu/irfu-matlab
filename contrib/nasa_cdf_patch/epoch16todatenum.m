function out = EPOCH16toDATENUM(epoch16)
%EPOCH16toDATENUM converts the time in UTC string (returned from cdfread) or
%                 date/time values for CDF_EPOCH16 to MATLAB datenum
%
%   OUT = EPOCH16toDATENUM(epoch16) returns MATLAB datenum.
%   OUT a column vector of numerical values of MATLAB date numbers.
%
%     epoch16               A vector or cell of UTC string or an M-by-10 matrix
%                           containing M full or partial date vectors for
%                           year, month, day, hour, minute, second, millisecond,
%                           microsecond, nanosecond and picosecond, in that 
%                           order.
%
%   Note:
%     The valid epoch string should have one of the following forms:
%     1. dd-mmm-yyyy hh:mm:ss.mmm.uuu.nnn.ppp (length of 36), e.g.,
%       "01-JAN-2000 12:00:00.123.456.789.000"
%     2. yyyymmdd.ddddddddddddddd (length of 24), e.g.,
%       "20000101.120000000000000"
%     3. yyyymmddhhmmss (length of 14), e.g.,
%       "20000101120000"
%     4. yyyy-mm-ddThh:mm:ss.mmmuuunnnppp (length of 32), e.g.,
%       "2000-01-01T12:00:00.123456789000"
%     where mmmuuunnnppp is milliseconds, microseconds, nanoseconds and
%                           picoseconds.
%
%   Examples:
%
%   % Read all the variable data in a CDF file. Among them, the variable of 
%     CDF_EPOCH16 data type is returned in UTC strings. Convert the strings
%     to MATLAB datenum. The variable is at index of 17 from the output of
%     1 by 20 array of cells from cdfread, the efficient way. 
%
%   data = cdfread('test', 'combinerecords', true);
%   epoch = data(1,17);
%   datenums = EPOCH16toDATENUM(epoch);
%
%   % Similar to the above example. But, read in the data into 40 by 20 array
%     of cells by cdfread. Convert the variable of CDF_EPOCH16 data type in
%     UTC strings and convert them to MATLAB datenum.
%
%   data = cdfread('test');
%   epoch = data(:,17);
%   datenums = EPOCH16toDATENUM(epoch);
%
%   % Convert the UTC strings in vector to MATLAB datenum.
%
%   epoch1 = ['2009-01-01T00:00:00.123456789000';
%             '2009-01-01T12:00:00.123456789000'];
%   datenums = EPOCH16toDATENUM(epoch1);
%
%   % Convert the date/times in matrix to MATLAB datenum.
%
%   epoch2 = [2009 01 01 00 00 00 123 456 789 000;
%             2009 01 01 12 00 00 123 456 789 000];
%   datenums = EPOCH16toDATENUM(epoch2);
%
%   See also CDFREAD.

%   Copyright 1984-2009 The MathWorks, Inc.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:EPOCH16toDATENUM:inputArgumentCount', ...
          'EPOCH16toDATENUM requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:EPOCH16toDATENUM:outputArguments', ...
          'EPOCH16toDATENUM requires only one output argument.')
end
s = size(epoch16);
if (iscell(epoch16(1,1)))
  if (s(1) == 1)
    epoch16 = epoch16{1,1};
  else
    if (iscellstr(epoch16))
      epoch16 = char(epoch16);
    else
      if (s(2) == 1)
        for x=1:s(1)
        epoch16(x,1) = epoch16{x,1};
        end
      end
    end
  end
end
if (iscell(epoch16))
  epoch16 = cell2mat(epoch16);
end
if (isnumeric(epoch16))
  s = size(epoch16);
  if (s(2) ~= 10)
    error('MATLAB:EPOCH16toDATENUM:inputArgumentvector', ...
          'EPOCH16toDATENUM requires an M by 10 vector for date/time fields.')
  end
  dates=epoch16(:,1:6);
  for p =1:s(1)
    dates(p,6) = epoch16(p,6) + epoch16(p,7)/1000 + epoch16(p,8)/1000000 + ...
                 epoch16(p,9)/1000000000 + epoch16(p,10)/1000000000000;
  end
  out = datenum(dates); 
else
  dates = parseepoch16c(epoch16);
  out = datenum(dates);
end
