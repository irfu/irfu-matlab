function out = SPDFEPOCH16UnixTime(time, varargin)
%SPDFEPOCH16UnixTime converts the CDF_EPOCH16 time (picoseconds since 
%                  01-01-0000 00:00:00.000) to unix time (seconds from 
%                  01-01-1970 00:00:00.000) or vice verse. The Unix time can
%                  have sub-second, with a time resolution of microseconds, 
%                  in its fractional part.
%
%   OUT = SPDFEPOCH16UnixTime(time, 'TOEPOCH16', TF) returns converted time(s).
%   OUT a column vector of numerical values of converted times.
%
%     time         A vector or scalar of time(s), based on CDF_EPOCH16 data 
%                  type or unix time. 
%
%   The option 'TOEPOCH16' is specified to true if the time conversion is from
%   unix time to CDF_EPOCH16. If false, or not specified, the conversion is
%   from CDF_EPOCH16 to unix time.
%
%   Note: Since CDF_EPOCH16 has higher time resolution, the converted
%         time might not preserve the full resolution for sub-microseconds.
%         sub-microseconds portion.
%
%   Examples:
%
%   % Convert a CDF_EPOCH16 data, a scalar at 
%     01-01-2000 00:00:00.123.456.789.000 to a unix
%     time value.
%
%   epoch = spdfcomputeepoch16([2000,1,1,0,0,0,123,456,789,000]);
%   unixtime = SPDFEPOCH16UnixTime(epoch);
%
%   % Convert CDF_EPOCH16 data, a vector of
%     01-01-2000 00:00:00.123.456.789.000 and 
%     02-02-2001 01:01:01.456.789.012.345 to unix time values.
%
%   epochs = spdfcomputeepoch16([2000,1,1,0,0,0,123,456,789,0;
%                                2001,2,2,1,1,1,456,789,012,345]);
%   unixtimes = SPDFEPOCH16UnixTime(epochs);
%
%   % Convert a unix time to CDF_EPOCH16 epoch. 
%
%   epoch = SPDFEPOCH16UnixTime(unixtime, 'TOEPOCH16', true);
%
%   % Convert a vector of unix times to CDF_EPOCH16 epochs.
%
%   epochs = SPDFEPOCH16UnixTime(unixtimes, 'TOEPOCH16', true);
%

% HISTORY:
%   November 22, 2018  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFEPOCH16UnixTime:inputArgumentCount', ...
          'SPDFEPOCH16UnixTime requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFEPOCH16UnixTime:outputArguments', ...
          'SPDFEPOCH16UnixTime requires only one output argument.')
end

[args, msg] = parse_inputs(varargin{:});

if (~isempty(msg))
    error('MATLAB:spdfepochunixtime:badInputArguments', '%s', msg)
end
if (args.toepoch == 0)
  time = transpose(time);
end
out = spdfepoch16unixtimec(time, args.toepoch);
if (args.toepoch == 1)
  out = transpose(out);
end

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
    paramStrings = {'toepoch16'};

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
         case 'toepoch16'
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
                   if (convert > 0)
                     args.toepoch = int32(1);
                   else
                     args.toepoch = int32(0);
                   end
               else
                   msg = 'Epoch conversion value must be a scalar logical.';
               end
           end
       end  % switch
    end  % for

end  % if (nargin > 1)

end
