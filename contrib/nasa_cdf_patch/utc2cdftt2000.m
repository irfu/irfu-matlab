function out = UTC2CDFTT2000(UTC,month,day,hour,minute,second,milsec,micsec,nansec)
%UTC2CDFTT2000 converts a UTC date/time in string or components to CDF_TIME_TT2000
%
%   There are two forms of this function:
%
%   OUT = UTC2CDFTT2000(UTC) 
%         Parses a single CDF epoch string and returns a CDF_TIME_TT2000
%         data type.
%
%     UTC               A UTC string
%
%   Note: the valid epoch string should be one of the following forms:
%      0: dd-mmm-yyyy hh:mm:ss.llluuunnn, e.g., "01-JAN-2010 12:00:00.000000000"
%      1: yyyymmdd.dddddddddd, e.g., "20100101.1200000000"
%      2: yyyymmddhhmmss, e.g., "20100101120000"
%      3: yyyy-mm-ddThh:mm:ss.llluuunnn, e.g., "2010-01-01T12:00:00.000000000"
%         where lll as milliseconds, uuu as microseconds and nnn as nanoseconds.
%   Or,
%
%   OUT = UTC2CDFTT2000(year,month,day,hour,minute,second,milsec,micsec,nansec) 
%         Compute the CDF epoch in CDF_TIME_TT2000 data type.
%
%     year, month, day, hour,       Integer form for date/time components
%     minute, second, milsec,  
%     micsec, nansec
%
%   Examples:
%
%   % Write 100 epoch records, starting from 2009-01-01T00:00:00.123456789 with
%   % one (1) second stride, to a CDF.
%
%   utc = '2009-01-01T00:00:00.123456789';
%   epoch = UTC2CDFTT2000(utc);
%   epochs = epoch+[0:99]*1000;
%   SPDFCDFWRITE('example', {'Epoch', epochs}, 'EPOCHArrayisCDFEpoch', true, ...
%            'recordbound', {'Epoch'});
%
%   % Alternatively, the sample can be entered as
%
%   epoch = UTC2CDFTT2000(2010,1,1,0,0,0,0,0,0);
%   epochs = epoch:[0:99]*1000;
%   SPDFCDFWRITE('example', {'Epoch', epochs}, 'EPOCHArrayisCDFEpoch', true, ...
%            'recordbound', {'Epoch'});
%
%   See also SPDFCDFWRITE, SPDFCDFREAD, CDFTT2000.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin ~= 1 && nargin ~= 9)
    error('MATLAB:UTC2CDFTT2000:inputArgumentCount', ...
          'UTC2CDFTT2000 requires one input string or nine date/time components.')
end

if (nargout > 1)
    error('MATLAB:UTC2CDFTT2000:outputArguments', ...
          'UTC2CDFTT2000 requires only one output argument.')
end

if (nargin == 1 && ischar(UTC))
  out = spdfparsett2000(UTC);
elseif (nargin == 9) 
  out = spdfcomputett2000([UTC,month,day,hour,minute,second,milsec,micsec,nansec]);
else
  error('MATLAB:UTC2CDFTT2000:input', ...
        'UTC2CDFTT2000 input does not meet requirements.')
end
end
