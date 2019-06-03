function out = SPDFPARSETT2000(tt2000)
%SPDFPARSETT2000 converts the CDF epoch in UTC string to CDF_TIME_TT2000
%data type.
%
%   OUT = SPDFPARSETT2000(tt2000) returns the CDF epoch in TT2000 data type.
%   OUT a vector of numerical values of mxINT64_CLASS (int64).
%
%     tt2000               A cell or vector of TT2000 UTC string
%
%   Note: the valid epoch string should be one of the following forms:
%      0: dd-mmm-yyyy hh:mm:ss.mmmuuunnn, e.g., "01-JAN-2010 12:00:00.000000000"
%      1: yyyymmdd.dddddddddd, e.g., "20100101.1200000000"
%      2: yyyymmddhhmmss, e.g., "20100101120000"
%      3: yyyy-mm-ddThh:mm:ss.mmmuuunnn, e.g., "2010-01-01T12:00:00.000000000"
%         where mmmuuunnn is for milliseconds, microseconds and nanoseconds.
%      4: yyyy-mm-ddThh:mm:ss.mmmuuunnnZ, e.g., "2010-01-01T12:00:00.000000000Z"
%         where mmmuuunnn is for milliseconds, microseconds and nanoseconds.
%
%   Examples:
%
%   % Convert the UTC strings in cell to TT2000 and write it to Epoch
%   % variable of CDF TT2000 data type in a CDF.
%
%   utcs = {'2009-01-01T00:00:00.123456789'; '2009-01-01T12:00:00.123456789'};
%   tt2000 = SPDFPARSETT2000(utcs);
%   SPDFCDFWRITE('example', {'Epoch', tt2000}, 'TT2000', true, ...
%            'recordbound', {'Epoch'});
%
%   See also CDFTT2000, SPDFENCODETT2000, SPDFCOMPUTETT2000, SPDFBREAKDOWNTT2000

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFPARSETT2000:inputArgumentCount', ...
          'SPDFPARSETT2000 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFPARSETT2000:outputArguments', ...
          'SPDFPARSETT2000 requires only one output argument.')
end

out = spdfparsett2000c(tt2000);

