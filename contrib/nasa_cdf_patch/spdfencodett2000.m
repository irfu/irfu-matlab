function out = SPDFENCODETT2000(tt2000, varargin)
%SPDFENCODETT2000 encodes the CDF epoch data in TT2000 data type to UTC
%strings.
%
%   OUT = SPDFENCODETT2000(tt2000) returns the CDF epoch in string. OUT is
%   a cell of strings.
%
%     tt2000               A single or vector of numerical values of
%                          CDF_TIME_TT2000 (an mxINT64_CLASS) data type.
%
%   The encoded epoch string will have the following format:
%       yyyy-mm-ddThh:mm:ss.mmmuuunnn, e.g., "2000-01-01T12:34:56.123456789"
%
%   OUT = SPDFENCODETT2000(tt2000, 'Format', FORMAT) encodes the UTC
%   string into the specified format. FORMAT is a number from 0 to 4.
%   FORMAT:
%     0: dd-mmm-yyyy hh:mm:ss.mmmuuunnn, e.g., "01-JAN-2000 12:34:56.123456789"
%     1: yyyymmdd.dddddddddd, e.g., "20000101.1200000000"
%     2: yyyymmddhhmmss, e.g., "20000101123456"
%     3: yyyy-mm-ddThh:mm:ss.mmmuuunnn, e.g., "2000-01-01T12:34:56.123456789"
%     4: yyyy-mm-ddThh:mm:ss.mmmuuunnnZ, e.g., "2000-01-01T12:34:56.123456789Z"
%   where mmmuuunnn is milliseconds, microseconds and nanoseconds.
%   Format 3 is the default (as ISO 8601) if this option is not provided.
%
%   Examples:
%
%   % Read all of the data from the same file, the most efficient way,
%   % and encode the Variable Epoch, of CDF_TIME_TT2000 data type, to UTC
%   % string.
%
%   data = spdfcdfread('example', 'Variables', {'Epoch'}, "KEEPEPOCHASIS", true);
%   out = SPDFENCODETT2000(data);
%
%
%   See also CDFTT2000, SPDFPARSETT2000, SPDFCOMPUTETT2000, SPDFBREAKDOWNTT2000,

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.
%   October 6, 2018  Mike Liu    Added format 4.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFENCODETT2000:inputArgumentCount', ...
          'SPDFENCODETT2000 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFENCODETT2000:outputArguments', ...
          'SPDFENCODETT2000 requires only one output argument.')
end

[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:SPDFENCODETT2000:badInputArguments', '%s', msg)
end
if (isa(tt2000,'cdftt2000'))
  if (length(tt2000) > 1)
    for p = 1:length(tt2000)
      if (~isa(tt2000, 'cell'))
        dataaa(p,1) = spdfencodett2000c(int64(todatenum(tt2000(p,1))), args.Format);
      else
        if (length(tt2000{p}) > 1)
          for q = 1:length(tt2000{p})
            dataaa(q) = spdfencodett2000c(int64(todatenum(tt2000{p}(q,1))), args.Format);
          end
        else
          dataaa(p,1) = spdfencodett2000c(int64(todatenum(tt2000{p})), args.Format);
        end
      end
    end
    out = dataaa;
  end
elseif (isa(tt2000,'double'))
  out = datestr(tt2000);
else
  out = spdfencodett2000c(tt2000, args.Format);
end
%%%
%%% Function parse_inputs
%%%

function [args, msg] = parse_inputs(varargin)
% Set default values
args.Format = int32(3);
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
                 msg = sprintf('format value "%d" is out of 0-4 range.', format);
                 return;
               end
               args.Format = int32(format);
           end
       end  % switch
    end  % for 
                   
end  % if (nargin > 1)


