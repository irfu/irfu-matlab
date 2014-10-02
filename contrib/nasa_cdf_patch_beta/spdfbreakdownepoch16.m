function out = SPDFBREAKDOWNEPOCH16(epoch16)
%SPDFBREAKDOWNEPOCH16 converts the CDF_EPOCH16 time, picoseconds since
%             0000-01-01 to UTC date/time.
%
%   OUT = SPDFBREAKDOWNEPOCH16(epoch16) returns the UTC date/time from CDF_EPOCH16
%   time. OUT is an array with each row having ten (10) numerical values
%   for year, month, day, hour, minute, second, millisecond, microsecond,
%   nanosecond and picosecond..
%
%     epoch16            An array with each row having two (2) numerical
%                        value in CDF_EPOCH16 of mxDOUBLE_CLASS (double).
%
%   Examples:
%
%   % breakdown three CDF_EPOCH16 times to their UTC form.
%
%   data = spdfparseepoch16 (['2009-01-01T00:00:00.123456789123', ...
%                             '2009-01-01T12:00:00.987654321123', ...
%                             '2009-01-01T12:34:56.123456789123']);
%   out = SPDFBREAKDOWNEPOCH16(data)
%   ans =
%
%        2009     1     1     0      0     0    123    456    789  123
%        2009     1     1    12      0     0    987    654    321  123
%        2010     1     1    12     34    56    123    456    789  123
%
%
%   See also SPDFENCODEEPOCH16, SPDFPARSEEPOCH16, SPDFCOMPUTEEPOCH16


% HISTORY:
%   August 16, 2014  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFBREAKDOWNEPOCH16:inputArgumentCount', ...
          'SPDFBREAKDOWNEPOCH16 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFBREAKDOWNEPOCH16:outputArguments', ...
          'SPDFBREAKDOWNEPOCH16 requires only one output argument.')
end

if (~isa(epoch16,'numeric'))
    error('MATLAB:SPDFBREAKDOWNEPOCH16:inputArguments', ...
          'SPDFBREAKDOWNEPOCH16 requires the input to be in numeric.')
end

ss=size(epoch16);
if (ss(1) < 1 || ss(2) ~= 2)
    error('MATLAB:SPDFBREAKDOWNEPOCH16:inputArguments', ...
          'SPDFBREAKDOWNEPOCH16 requires each row to have just two (2) field.')
end

  out = spdfbreakdownepoch16c (epoch16);
end
