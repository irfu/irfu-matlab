function out = SPDFCOMPUTEEPOCH(datetime)
%SPDFCOMPUTEEPOCH converts the UTC date/time components to CDF_EPOCH time,
%              in milliseconds since 00-01-01.
%
%   OUT = SPDFCOMPUTEEPOCH(datetime) returns the CDF_EPOCH time. OUT is
%   a vector of double values of mxDOUBLE_CLASS (double).
%
%     datetime             An array with each row having seven (7) numerical
%                          values for year, month, day, hour, minute, second,
%                          millisecond.
%
%   Examples:
%
%   % Compute two UTC times into CDF_EPOCH.
%
%   data = [2009 01 01  0 0 0 123;
%           2009 01 01 12 0 0 987];
%   out = SPDFCOMPUTEEPOCH(data);
%
%  SPDFENCODEEPOCH(out) will show:
%   ans = 
%
%       '01-JAN-2009 00:00:00.123'
%       '01-JAB-2009 12:00:00.987'
%
%
%   See also CDFEPOCH, SPDFENCODEEPOCH, SPDFPARSEEPOCH, SPDFBREAKDOWNEPOCH

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFCOMPUTEEPOCH:inputArgumentCount', ...
          'SPDFCOMPUTEEPOCH requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFCOMPUTEEPOCH:outputArguments', ...
          'SPDFCOMPUTEEPOCH requires only one output argument.')
end

if (~isa(datetime,'numeric'))
    error('MATLAB:SPDFCOMPUTEEPOCH:inputArguments', ...
          'SPDFCOMPUTEEPOCH requires the input to be in numeric.')
end

ss=size(datetime);
if (ss(1) < 1 || ss(2) ~= 7)
    error('MATLAB:SPDFCOMPUTEEPOCH:inputArguments', ...
          'SPDFCOMPUTEEPOCH requires each row to have seven (7) date/time fields.')
end

out = spdfcomputeepochc (datetime);

