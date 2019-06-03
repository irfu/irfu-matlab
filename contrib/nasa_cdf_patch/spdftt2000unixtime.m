function out = SPDFTT2000UnixTime(time, varargin)
%SPDFTT2000UnixTime converts the CDF_TIME_TT2000 time (nanoseconds since 
%                  2000-01-01 12:00:00.000 with leap seconds) to unix time
%                  (seconds from 1970-01-01 00:00:00.000) or vice verse.
%                  The Unix time can have sub-second, with a time resolution
%                  of microseconds, in its fractional part.
%
%   OUT = SPDFTT2000UnixTime(time, 'TOTT2000', TF) returns converted time(s).
%   OUT a column vector of numerical values of converted times.
%
%     time         A vector or scalar of time(s), based on CDF_TIME_TT2000 data type
%                  or unix time. 
%
%   The option 'TOTT2000' is specified to true if the time conversion is from
%   unix time to CDF_TIME_TT2000. If false, or not specified, the conversion is
%   from CDF_TIME_TT2000 to unix time.
%
%   Note: Since unix time does not include leap seconds, the time conversion will
%         not be properly handled if a time falls on a leap second. The TT2000
%         time has a higher time resolution, sub-microseonds might not be preserved
%         in unix time.
%
%   Examples:
%
%   % Convert a CDF_TIME_TT2000 data, a scalar at 2000-01-01 00:00:00.123 to a unix
%     time value.
%
%   epoch = spdfcomputett2000([2000,1,1,0,0,0,123,0,0]);
%   unixtime = SPDFTT2000UnixTime(epoch);
%
%   % Convert CDF_TIME_TT2000 data, a vector at 2000-01-01-2000 00:00:00.123 and 
%     2001-02-02 01:01:01.456 to unix time values.
%
%   epochs = spdfcomputett2000([2000,1,1,0,0,0,123,0,0; 2001,2,2,1,1,1,456,0,0]);
%   unixtimes = SPDFTT2000UnixTime(epochs);
%
%   % Convert a unix time to CDF_TIME_TT2000 epoch. 
%
%   epoch = SPDFTT2000UnixTime(unixtime, 'TOTT2000', true);
%
%   % Convert a vector of unix times to CDF_TIME_TT2000 epochs.
%
%   epochs = SPDFTT2000UnixTime(unixtimes, 'TOTT2000', true);
%

% HISTORY:
%   November 22, 2018  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFTT2000UnixTime:inputArgumentCount', ...
          'SPDFTT2000UnixTime requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFTT2000UnixTime:outputArguments', ...
          'SPDFTT2000UnixTime requires only one output argument.')
end

[args, msg] = parse_inputs(varargin{:});

if (~isempty(msg))
    error('MATLAB:spdftt2000unixtime:badInputArguments', '%s', msg)
end

out = spdftt2000unixtimec(time, args.toepoch);

end

%%%
%%% Function parse_inputs
%%%

function [args, msg] = parse_inputs(varargin)
% Set default values
args.toepoch = int32(0);
msg = '';
% Parse arguments based on their number.
if (nargin > 0)
    paramStrings = {'tott2000'};

    % For each pair
    for k = 1:2:length(varargin)
       param = lower(varargin{k});
       idx = strmatch(param, paramStrings);
       if (isempty(idx))
           msg = sprintf('Unrecognized parameter name "%s".', param);
           return
       elseif (length(idx) > 1)
           msg = sprintf('Ambiguous parameter name "%s".', param);
           return
       end
       switch (paramStrings{idx})
         case 'tott2000'
           if (k == length(varargin))
               msg = 'No epoch conversion value specified.';
               return
           else
               convert = varargin{k + 1};
               if (numel(convert) ~= 1)
                   msg = 'Conversion value must be a scalar logical.';
               end
               if (islogical(convert))
                   if (convert)
                     args.toepoch = int32(1);
                   else
                     args.toepoch = int32(0);
                   end
               elseif (isnumeric(convert))
                   args.toepoch = int32(convert);
               else
                   msg = 'Epoch conversion value must be a scalar logical.';
               end
           end
       end  % switch
    end  % for

end  % if (nargin > 1)

end
