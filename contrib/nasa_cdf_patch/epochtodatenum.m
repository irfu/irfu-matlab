function out = EPOCHtoDATENUM(epoch)
%EPOCHtoDATENUM converts the time in UTC string (returned from cdfread) or
%               date/time values for CDF_EPOCH to MATLAB datenum
%
%   OUT = EPOCHtoDATENUM(epoch) returns MATLAB datenum.
%   OUT a column vector of numerical values of MATLAB date numbers.
%
%     epoch               A vector or cell of UTC string or an M-by-7 matrix
%                           containing M full or partial date vectors for
%                           year, month, day, hour, minute, second, millisecond,
%                           in that order.
%
%   Note:
%     The valid epoch string should have one of the following forms:
%     1. dd-mmm-yyyy hh:mm:ss.lll (length of 24), e.g.,
%       "01-JAN-2000 12:00:00.123"
%     2. yyyymmdd.ddddddd (length of 16), e.g.,
%       "20000101.1200000"
%     3. yyyymmddhhmmss (length of 14), e.g.,
%       "20000101120000"
%     4. yyyy-mm-ddThh:mm:ss.lll (length of 23), e.g.,
%       "2000-01-01T12:00:00.123"
%
%   Examples:
%
%   % Read all the variable data in a CDF file. Among them, the variable of 
%     CDF_EPOCH data type is returned in UTC strings. Convert the strings
%     to MATLAB datenum. The variable is at index of 17 from the output of
%     1 by 20 array of cells from cdfread, the efficient way. 
%
%   data = cdfread('test', 'CombineRecords', true, 'KeepEpochAsIs', true);
%   epoch = data(1,17);
%   datenums = EPOCHtoDATENUM(epoch);
%
%   % Similar to the above example. But, read in the data into 40 by 20 array
%     of cells by cdfread. Convert the variable of CDF_EPOCH data type in
%     UTC strings and convert them to MATLAB datenum.
%
%   data = cdfread('test');
%   epoch = data(:,17);
%   datenums = EPOCHtoDATENUM(epoch);
%
%   % Convert the UTC strings in vector to MATLAB datenum.
%
%   epoch1 = ['2009-01-01T00:00:00.123';
%             '2009-01-01T12:00:00.123'];
%   datenums = EPOCHtoDATENUM(epoch1);
%
%   % Convert the date/times in matrix to MATLAB datenum.
%
%   epoch2 = [2009 01 01 00 00 00 123;
%             2009 01 01 12 00 00 123];
%   datenums = EPOCHtoDATENUM(epoch2);
%
%   See also CDFREAD.

%   Copyright 1984-2009 The MathWorks, Inc.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:EPOCHtoDATENUM:inputArgumentCount', ...
          'EPOCHtoDATENUM requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:EPOCHtoDATENUM:outputArguments', ...
          'EPOCHtoDATENUM requires only one output argument.')
end
s = size(epoch);
if (iscell(epoch(1,1)))
  if (s(1) == 1)
    epoch = epoch{1,1};
  else
    if (iscellstr(epoch))
      epoch = char(epoch);
    else
      if (s(2) == 1)
        for x=1:s(1)
        epoch(x,1) = epoch{x,1};
        end
      end
    end
  end
end
if (iscell(epoch))
  epoch = cell2mat(epoch);
end
if (isnumeric(epoch))
  s = size(epoch);
  if (s(2) == 1)
    out = datenum(epoch./86400000.0+1);
  else
    if (s(2) ~= 7)
      error('MATLAB:EPOCHtoDATENUM:inputArgumentvector', ...
            'EPOCHtoDATENUM requires an M by 7 vector for date/time fields.')
    end
    dates=epoch(:,1:6);
    for p =1:s(1)
      dates(p,6) = epoch(p,6) + epoch(p,7)/1000.0;
    end
    out = datenum(dates); 
  end
else
  dates = parseepochc(epoch);
  out = datenum(dates);
end
