function out = ENCODETT2000(tt2000, varargin)
%ENCODETT2000 encodes the CDF epoch data in TT2000 data type to UTC
%strings.
%
%   OUT = ENCODETT2000(tt2000) returns the CDF epoch in string. OUT is
%   a cell of strings.
%
%     tt2000               A single or vector of numerical values of
%                          CDF_TIME_TT2000 (an mxINT64_CLASS) data type.
%
%   The encoded epoch string will have the following format:
%       yyyy-mm-ddThh:mm:ss.mmmuuunnn, e.g., "2000-01-01T12:34:56.123456789"
%
%   OUT = ENCODETT2000(tt2000, 'Format', FORMAT) encodes the UTC
%   string into the specified format. FORMAT is a number from 0 to 3.
%   FORMAT:
%     0: dd-mmm-yyyy hh:mm:ss.mmmuuunnn, e.g., "01-JAN-2000 12:34:56.123456789"
%     1: yyyymmdd.dddddddddd, e.g., "20000101.1200000000"
%     2: yyymmddhhmmss, e.g., "20000101123456"
%     3: yyyy-mm-ddThh:mm:ss.mmmuuunnn, e.g., "2000-01-01T12:34:56.123456789"
%   where mmmuuunnn is milliseconds, microseconds and nanoseconds.
%   Format 3 is the default (as ISO 8601) if this option is not provided.
%
%   Examples:
%
%   % Read all of the data from the file, example.cdf, the inefficient way,
%   % and encode the TT2000 data type to UTC string at column 19.
%
%   data = cdfread('example');
%   out = ENCODETT2000(data(:,19));
%
%   % Read all of the data from the same file, the most efficient way,
%   % and encode the Variable Epoch, of CDF_TIME_TT2000 data type, to UTC
%   % string.
%
%   data = cdfread('example', 'CombineRecords', true, 'Variables', {'Epoch'});
%   out = ENCODETT2000(data);
%
%
%   See also CDFTT2000, CDFREAD, PARSETT2000, COMPUTETT2000, BREAKDOWNTT2000,
%            CDFLEAPSECONDSINFO, ENCODEEPOCH.

%   Copyright 1984-2009 The MathWorks, Inc.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:ENCODETT2000:inputArgumentCount', ...
          'ENCODETT2000 requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:ENCODETT2000:outputArguments', ...
          'ENCODETT2000 requires only one output argument.')
end

[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:ENCODETT2000:badInputArguments', '%s', msg)
end
if (isa(tt2000,'cdftt2000'))
  if (length(tt2000) > 1)
    for p = 1:length(tt2000)
      if (~isa(tt2000, 'cell'))
        dataaa(p,1) = encodett2000c(int64(todatenum(tt2000(p,1))), args.Format);
      else
        if (length(tt2000{p}) > 1)
          for q = 1:length(tt2000{p})
            dataaa(q) = encodett2000c(int64(todatenum(tt2000{p}(q,1))), args.Format);
          end
        else
          dataaa(p,1) = encodett2000c(int64(todatenum(tt2000{p})), args.Format);
        end
      end
    end
    out = dataaa;
  end
elseif (isa(tt2000,'double'))
  out = datestr(tt2000);
else
  out = encodett2000c(tt2000, args.Format);
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

               if (int32(format) < 0 || int32(format) > 3)
                 msg = sprintf('format value "%d" is out or 0-3 range.', format);
                 return;
               end
               args.Format = int32(format);
           end
       end  % switch
    end  % for 
                   
end  % if (nargin > 1)


