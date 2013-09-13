function out = BREAKDOWNTT2000(tt2000)
%BREAKDOWNTT2000 converts the CDF TT2000 time, nanoseconds since
%             2000-01-01 12:00:00 to UTC date/time.
%
%   OUT = BREAKDOWNTT2000(tt2000) returns the UTC date/time from CDF TT2000
%   time. OUT is an array with each row having nine (9) numerical values
%   for year, month, day, hour, minute, second, millisecond, microsecond
%   and nanosecond.
%
%     tt2000             An array with each row having one (1) numerical
%                        value in CDF TT2000 of mxINT64_CLASS (int64).
%
%   Examples:
%
%   % breakdown three CDF TT2000 times to their UTC form.
%
%   data = [ 284040066307456789; 284083267171654321; 315621362307456789];
%   out = BREAKDOWNTT2000(data)
%   ans =
%
%        2009     1     1     0      0     0    123    456    789
%        2009     1     1    12      0     0    987    654    321
%        2010     1     1    12     34    56    123    456    789
%
%
%   See also CDFTT2000, CDFREAD, ENCODETT2000, PARSETT2000, COMPUTETT2000,
%            CDFLEAPSECONDSINFO.

%   Copyright 1984-2009 The MathWorks, Inc.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:BREAKDOWNTT2000:inputArgumentCount', ...
          'BREAKDOWNTT2000 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:BREAKDOWNTT2000:outputArguments', ...
          'BREAKDOWNTT2000 requires only one output argument.')
end

if (~isa(tt2000,'numeric'))
    error('MATLAB:BREAKDOWNTT2000:inputArguments', ...
          'BREAKDOWNTT2000 requires the input to be in numeric.')
end

ss=size(tt2000);
if (ss(1) < 1 || ss(2) ~= 1)
    error('MATLAB:BREAKDOWNTT2000:inputArguments', ...
          'BREAKDOWNTT2000 requires each row to have just one (1) field.')
end

if (~isa(tt2000, 'int64'))
  out = breakdowntt2000c (int64(tt2000(:)));
else
  out = breakdowntt2000c (tt2000);
end
