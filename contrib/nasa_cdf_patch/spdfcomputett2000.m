function out = SPDFCOMPUTETT2000(datetime)
%SPDFCOMPUTETT2000 converts the UTC date/time components to CDF TT2000 time,
%              in nanoseconds since 2000-01-01 12:00:00 with leap seconds.
%
%   OUT = SPDFCOMPUTETT2000(datetime) returns the CDF TT2000 time. OUT is
%   a vector of integer values of mxINT64_CLASS (int64).
%
%     datetime             An array with each row having nine (9) numerical
%                          values for year, month, day, hour, minute, second,
%                          millisecond, microsecond and nanosecond.
%
%   Examples:
%
%   % Compute two UTC times into CDF TT2000.
%
%   data = [2009 01 01  0 0 0 123 456 789;
%           2009 01 01 12 0 0 987 654 321];
%   out = SPDFCOMPUTETT2000(data);
%
%   SPDFENCODETT2000(out) will show:
%   ans = 
%
%       '2009-01-01T00:00:00.123456789'
%       '2009-01-01T12:00:00.987654321'
%
%
%   See also CDFTT2000, SPDFENCODETT2000, SPDFPARSETT2000, SPDFBREAKDOWNTT2000

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFCOMPUTETT2000:inputArgumentCount', ...
          'SPDFCOMPUTETT2000 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFCOMPUTETT2000:outputArguments', ...
          'SPDFCOMPUTETT2000 requires only one output argument.')
end

if (~isa(datetime,'numeric'))
    error('MATLAB:SPDFCOMPUTETT2000:inputArguments', ...
          'SPDFCOMPUTETT2000 requires the input to be in numeric.')
end

ss=size(datetime);
if (ss(1) < 1 || ss(2) ~= 9)
    error('MATLAB:SPDFCOMPUTETT2000:inputArguments', ...
          'SPDFCOMPUTETT2000 requires each row to have nine (9) date/time fields.')
end

out = spdfcomputett2000c (datetime);

