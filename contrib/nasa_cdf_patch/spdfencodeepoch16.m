function out = SPDFENCODEEPOCH16(epoch, varargin)
%SPDFENCODEEPOCH16 encodes an epoch of CDF_EPOCH16 data type, an double array value
%
%   OUT = SPDFENCODEEPOCH16(epoch) 
%         Returns a UTC string(s).
%
%     epoch                An epoch
%
%   The encoded epoch string will have the following ISO 8601 format:
%          yyyy-mm-ddThh:mm:ss.mmmuuunnnppp
%          e.g., "2000-01-01T12:34:56.123456789999"
%   Originally, it was in this form:
%          dd-mmm-yyyy hh:mm:ss.mmm.uuu.nnn.ppp
%          e.g., "01-Jan-2000 12:34:56.123.456.789.999"
%
%   OUT = SPDFENCODEEPOCH16(epoch, 'Format', FORMAT) encodes the UTC
%   string into the specified format. FORMAT is a number from 0 to 4.
%   FORMAT:
%     0: dd-mmm-yyyy hh:mm:ss.mmm.uuu.nnn.ppp 
%        e.g., "01-JAN-2000 12:34:56.123.456.789.999"
%     1: yyyymmdd.ddddddd e.g., "20000101.1200000"
%     2: yyyymmddhhmmss e.g., "20000101123456"
%     3: yyyy-mm-ddThh:mm:ss.mmm.uuu.nnn.pppZ
%        e.g., "2000-01-01T12:34:56.123.456.789.999Z"
%     4: yyyy-mm-ddThh:mm:ss.mmmuuunnnppp
%        e.g., "2000-01-01T12:34:56.123456789999"
%   where mmm is milliseconds,
%   where uuu is microseconds,
%   where nnn is nanoseconds,
%   where ppp is picoseconds.
%
%   Examples:
%
%   % Acquire 'Epoch' variable data as is (two double values) from 'sample' CDF
%   % and encode the epoch values.
%
%   epochs = spdfcdfread('Sample', 'Variables', {'Epoch'}, 
%                        "KeepEpochAsIs", true);
%   spdfencodeepoch16(epochs)
%
%   See also SPDFBREAKDOWNEPOCH16, SPDFCOMPUTEEPOCH16, SPDFPARSEEPOCH16

% HISTORY:
%   March 5, 2013  Mike Liu    The initial version.

if (nargin < 1)
    error('MATLAB:SPDFENCODEEPOCH16:inputArgumentCount', ...
          'SPDFENCODEEPOCH16 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFENCODEEPOCH16:outputArguments', ...
          'SPDFENCODEEPOCH16 requires only one output argument.')
end

[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:SPDFENCODEEPOCH:badInputArguments', '%s', msg)
end
out = spdfencodeepoch16c(epoch, args.Format);
%%%
%%% Function parse_inputs
%%%

function [args, msg] = parse_inputs(varargin)
% Set default values
args.Format = int32(4);
msg = '';
% Parse arguments based on their number.
if (nargin > 0)
    paramStrings = {'format'};

    % For each pair
    for k = 1:2:length(varargin)
       param = lower(varargin{k});
       if (~ischar(param))
           msg = 'Parameter name must be a string.';
           return
       end

       idx = strmatch(param, paramStrings);

       if (isempty(idx))
           msg = sprintf('Unrecognized parameter name "%s".', param);
           return
       elseif (length(idx) > 1)
           msg = sprintf('Ambiguous parameter name "%s".', param);
           return
       end

       switch (paramStrings{idx})
       case 'format'

           if (k == length(varargin))
               msg = 'No format specified.';
               return
           else
               format = varargin{k + 1};

               if (~isa(format, 'double') || ~isscalar(format))
                   msg = 'Format must be a single integer.';
                   return;
               end

               if (int32(format) < 0 || int32(format) > 4)
                 msg = sprintf('format value "%d" is out or 0-4 range.', format);
                 return;

               end
               args.Format = int32(format);
           end
       end  % switch
    end  % for

end  % if (nargin > 1)

