function out = SPDFPARSEEPOCH(epoch)
%SPDFPARSEEPOCH converts the epoch in UTC string to CDF_EPOCH
%data type.
%
%   OUT = SPDFPARSEEPOCH(epoch) returns the CDF epoch in EPOCH type.
%   OUT a vector of numerical values of mxDOUBLE_CLASS (double).
%
%     epoch               A cell or vector of UTC string
%
%   Note: the valid epoch string must be one of the following forms:
%      0: dd-mmm-yyyy hh:mm:ss.mmm, e.g., "01-JAN-2010 12:00:00.000"
%      1: yyyymmdd.ddddddd, e.g., "20100101.1200000"
%      2: yyyymmddhhmmss, e.g., "20100101120000"
%      3: yyyy-mm-ddThh:mm:ss.mmmZ, e.g., "2010-01-01T12:00:00.000Z"
%         where mmm is milliseconds.
%      4: yyyy-mm-ddThh:mm:ss.mmm, e.g., "2010-01-01T12:00:00.000"
%         where mmm is milliseconds.
%
%   Examples:
%
%   % Convert the UTC strings in cell to CDF_EPOCH and write it to Epoch
%   % variable of CDF CDF_EPOCH data type in a CDF.
%
%   utcs = {'2009-01-01 00:00:00.123'; '2009-01-01 12:00:00.123'};
%   epoch = SPDFPARSEEPOCH(utcs);
%   SPDFCDFWRITE('example', {'Epoch', epoch}, ...
%            'recordbound', {'Epoch'});
%
%   See also CDFEPOCH, SPDFENCODEEPOCH, SPDFCOMPUTEEPOCH, SPDFBREAKDOWNEPOCH

% HISTORY:
%   August 16, 2014  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFPARSEEPOCH:inputArgumentCount', ...
          'SPDFPARSEEPOCH requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFPARSEEPOCH:outputArguments', ...
          'SPDFPARSEEPOCH requires only one output argument.')
end

out = spdfparseepochc(epoch);

