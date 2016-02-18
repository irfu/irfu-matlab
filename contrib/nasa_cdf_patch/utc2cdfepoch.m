function out = UTC2CDFEPOCH(UTC,month,day,hour,minute,second,milsec,micsec)
%UTC2CDFEPOCH converts a UTC date/time in string or components to CDF_EPOCH
%
%   There are two forms of this function:
%
%   OUT = UTC2CDFEPOCH(UTC) 
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
%         where lll as milliseconds.
%   Or,
%
%   OUT = UTC2CDFEPOCH(year,month,day,hour,minute,second,milsec) 
%         Compute the CDF epoch in CDF_EPOCH data type.
%
%     year, month, day, hour,       Integer form for date/time components
%     minute, second, milsec.  
%
%   Examples:
%
%   % Write 100 epoch records, starting from 2009-01-01T00:00:00.123456789 with
%   % one (1) second stride, to a CDF.
%
%   utc = '2009-01-01T00:00:00.123';
%   epoch = UTC2CDFEPOCH(utc);
%   epochs = epoch+[0:99]*1000;
%   SPDFCDFWRITE('example', {'Epoch', epochs}, 'EPOCHArrayisCDFEpoch', true, ...
%            'recordbound', {'Epoch'});
%
%   % Alternatively, the sample can be entered as
%
%   epoch = UTC2CDFEPOCH(2010,1,1,0,0,0,0);
%   epochs = epoch:[0:99]*1000;
%   SPDFCDFWRITE('example', {'Epoch', epochs}, 'EPOCHArrayisCDFEpoch', true, ...
%            'recordbound', {'Epoch'});
%
%   See also SPDFCDFWRITE, SPDFCDFREAD, CDFEPOCH.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin ~= 1 && nargin ~= 7)
    error('MATLAB:UTC2CDFEPOCH:inputArgumentCount', ...
          'UTC2CDFEPOCH requires one input string or seven date/time components.')
end

if (nargout > 1)
    error('MATLAB:UTC2CDFEPOCH:outputArguments', ...
          'UTC2CDFEPOCH requires only one output argument.')
end

if (nargin == 1 && ischar(UTC))
  out = spdfparseepoch(UTC);
elseif (nargin == 7) 
  out = spdfcomputeepoch([UTC,month,day,hour,minute,second,milsec]);
else
  error('MATLAB:UTC2CDFEPOCH:input', ...
        'UTC2CDFEPOCH input does not meet requirements.')
end
end
