function out = SPDFBREAKDOWNEPOCH(epoch)
%SPDFBREAKDOWNEPOCH converts the CDF_EPOCH time, milliseconds since
%             0000-01-01 to UTC date/time.
%
%   OUT = SPDFBREAKDOWNEPOCH(epoch) returns the UTC date/time from CDF_EPOCH
%   time. OUT is an array with each row having seven (7) numerical values
%   for year, month, day, hour, minute, second, millisecond.
%
%     epoch            A vector with each row having one (1) numerical
%                      value in CDF_EPOCH of mxDOUBLE_CLASS (double).
%
%   Examples:
%
%   % breakdown three CDF_EPOCH times to their UTC form.
%
%   data = spdfparseepoch (['2009-01-01T00:00:00.123', ...
%                           '2009-01-01T12:00:00.987', ...
%                           '2009-01-01T12:34:56.123']);
%   out = SPDFBREAKDOWNEPOCH(data)
%   ans =
%
%        2009     1     1     0      0     0    123
%        2009     1     1    12      0     0    987
%        2010     1     1    12     34    56    123
%
%
%   See also CDFEPOCH, SPDFENCODEEPOCH, SPDFPARSEEPOCH, SPDFCOMPUTEEPOCH


% HISTORY:
%   August 16, 2014  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFBREAKDOWNEPOCH:inputArgumentCount', ...
          'SPDFBREAKDOWNEPOCH requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFBREAKDOWNEPOCH:outputArguments', ...
          'SPDFBREAKDOWNEPOCH requires only one output argument.')
end

ss=size(epoch);
if (ss(1) < 1 || ss(2) ~= 1)
    error('MATLAB:SPDFBREAKDOWNEPOCH:inputArguments', ...
          'SPDFBREAKDOWNEPOCH requires each row to have just one (1) field.')
end

  out = spdfbreakdownepochc (epoch);
end
