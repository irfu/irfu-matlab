function out = SPDFPARSEEPOCH16(epoch16)
%SPDFPARSEEPOCH16 converts the epoch in UTC string to CDF_EPOCH16
%data type.
%
%   OUT = SPDFPARSEEPOCH16(epoch16) returns the CDF epoch in CDF_EPOCH16 type.
%   OUT an array of numerical values of mxDOUBLE_CLASS (double).
%
%     epoch16               A cell or vector of UTC string
%
%   Note: the valid epoch string must be one of the following forms:
%      0: dd-mmm-yyyy hh:mm:ss.mmm.uuu.nnn.ppp, e.g., "01-JAN-2010 12:00:00.000.000.000.000"
%      1: yyyymmdd.ddddddddddddddd, e.g., "20100101.120000000000000"
%      2: yyyymmddhhmmss, e.g., "20100101120000"
%      3: yyyy-mm-ddThh:mm:ss.mmm.uuu.nnn.pppZ, e.g., "2010-01-01T12:00:00.000.000.000.000Z"
%         where mmm is milliseconds, uuu microseconds, nnn nanoseconds,
%         ppp picoseconds.
%      4: yyyy-mm-ddThh:mm:ss.mmmuuunnnppp, e.g., "2010-01-01T12:00:00.000000000000"
%         where mmm is milliseconds, uuu microseconds, nnn nanoseconds,
%         ppp picoseconds.
%
%   Examples:
%
%   % Convert the UTC strings in cell to CDF_EPOCH16 and write it to Epoch
%   % variable of CDF CDF_EPOCH16 data type in a CDF.
%
%   utcs = {'2009-01-01T00:00:00.123'; '2009-01-01T12:00:00.123'};
%   epoch16 = SPDFPARSEEPOCH16(utcs);
%   SPDFCDFWRITE('example', {'Epoch', epoch16}, ...
%            'recordbound', {'Epoch'});
%
%   See also SPDFENCODEEPOCH16, SPDFCOMPUTEEPOCH16, SPDFBREAKDOWNEPOCH16

% HISTORY:
%   August 16, 2014  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFPARSEEPOCH16:inputArgumentCount', ...
          'SPDFPARSEEPOCH16 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFPARSEEPOCH16:outputArguments', ...
          'SPDFPARSEEPOCH16 requires only one output argument.')
end

out = spdfparseepoch16c(epoch16);

