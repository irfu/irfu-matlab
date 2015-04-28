function out = SPDFBREAKDOWNTT2000(tt2000)
%SPDFBREAKDOWNTT2000 converts the CDF TT2000 time, nanoseconds since
%             2000-01-01 12:00:00 to UTC date/time.
%
%   OUT = SPDFBREAKDOWNTT2000(tt2000) returns the UTC date/time from CDF TT2000
%   time. OUT is an array with each row having nine (9) numerical values
%   for year, month, day, hour, minute, second, millisecond, microsecond
%   and nanosecond.
%
%     tt2000             A vector with each row having one (1) numerical
%                        value in CDF TT2000 of mxINT64_CLASS (int64).
%
%   Examples:
%
%   % breakdown three CDF TT2000 times to their UTC form.
%
%   data = [ 284040066307456789; 284083267171654321; 315621362307456789];
%   out = SPDFBREAKDOWNTT2000(data)
%   ans =
%
%        2009     1     1     0      0     0    123    456    789
%        2009     1     1    12      0     0    987    654    321
%        2010     1     1    12     34    56    123    456    789
%
%
%   See also CDFTT2000, SPDFENCODETT2000, SPDFPARSETT2000, SPDFCOMPUTETT2000

% HISTORY:
%   August 16, 2014  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFBREAKDOWNTT2000:inputArgumentCount', ...
          'SPDFBREAKDOWNTT2000 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFBREAKDOWNTT2000:outputArguments', ...
          'SPDFBREAKDOWNTT2000 requires only one output argument.')
end

if (~isa(tt2000,'numeric'))
    error('MATLAB:SPDFBREAKDOWNTT2000:inputArguments', ...
          'SPDFBREAKDOWNTT2000 requires the input to be in numeric.')
end

ss=size(tt2000);
if (ss(1) < 1 || ss(2) ~= 1)
    error('MATLAB:SPDFBREAKDOWNTT2000:inputArguments', ...
          'SPDFBREAKDOWNTT2000 requires each row to have just one (1) field.')
end

if (~isa(tt2000, 'int64'))
  out = spdfbreakdowntt2000c (int64(tt2000(:)));
else
  out = spdfbreakdowntt2000c (tt2000);
end
