function out = SPDFCOMPUTEEPOCH16(datetime)
%SPDFCOMPUTEEPOCH16 converts the UTC date/time components to CDF_EPOCH16 time,
%              in nanoseconds since 2000-01-01 12:00:00 with leap seconds.
%
%   OUT = SPDFCOMPUTEEPOCH16(datetime) returns the CDF_EPOCH16 time. OUT is
%   an array of N by 2 double values of mxDOUBLE_CLASS (double), each row
%   having two double values..
%
%     datetime             An array with each row having ten (10) numerical
%                          values for year, month, day, hour, minute, second,
%                          millisecond, microsecond and nanosecond.
%
%   Examples:
%
%   % Compute two UTC times into CDF_EPOCH16.
%
%   data = [2009 01 01  0 0 0 123 456 789 123;
%           2009 01 01 12 0 0 987 654 321 987];
%   out = SPDFCOMPUTEEPOCH16(data);
%
%   SPDFENCODEEPOCH16(out) will show:
%   ans = 
%
%       '01-JAN-2009 00:00:00.123456789123'
%       '01-JAN-2009 12:00:00.987654321987'
%
%
%   See also SPDFENCODEEPOCH16, SPDFPARSEEPOCH16, SPDFBREAKDOWNEPOCH16

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFCOMPUTEEPOCH16:inputArgumentCount', ...
          'SPDFCOMPUTEEPOCH16 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFCOMPUTEEPOCH16:outputArguments', ...
          'SPDFCOMPUTEEPOCH16 requires only one output argument.')
end

if (~isa(datetime,'numeric'))
    error('MATLAB:SPDFCOMPUTEEPOCH16:inputArguments', ...
          'SPDFCOMPUTEEPOCH16 requires the input to be in numeric.')
end

ss=size(datetime);
if (ss(1) < 1 || ss(2) ~= 10)
    error('MATLAB:SPDFCOMPUTEEPOCH16:inputArguments', ...
          'SPDFCOMPUTEEPOCH16 requires each row to have ten (10) date/time fields.')
end

out = spdfcomputeepoch16c (datetime);

