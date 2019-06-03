function out = SPDFEPOCHUnixTime(time, varargin)
%SPDFEPOCHUnixTime converts the CDF_EPOCH time (milliseconds since 
%                  01-01-0000 00:00:00.000) to unix time (seconds from 
%                  01-01-1970 00:00:00.000) or vice verse. The Unix time
%                  can have sub-second, with a time resolution of microseconds,
%                  in its fractional part.
%
%   OUT = SPDFEPOCHUnixTime(time, 'TOEPOCH', TF) returns converted time(s).
%   OUT a column vector of numerical values of converted times.
%
%     time         A vector or scalar of time(s), based on CDF_EPOCH data type
%                  or unix time. 
%
%   The option 'TOEPOCH' is specified to true if the time conversion is from
%   unix time to CDF_EPOCH. If false, or not specified, the conversion is
%   from CDF_EPOCH to unix time.
%
%   Examples:
%
%   % Convert a CDF_EPOCH data, a scalar at 01-01-2000 00:00:00.123 to a unix
%     time value.
%
%   epoch = spdfcomputeepoch([2000,1,1,0,0,0,123]);
%   unixtime = SPDFEPOCHUnixTime(epoch);
%
%   % Convert CDF_EPOCH data, a vector at 01-01-2000 00:00:00.123 and 
%     02-02-2001 01:01:01.456 to unix time values.
%
%   epochs = spdfcomputeepoch([2000,1,1,0,0,0,123; 2001,2,2,1,1,1,456]);
%   unixtimes = SPDFEPOCHUnixTime(epochs);
%
%   % Convert a unix time to CDF_EPOCH epoch. 
%
%   epoch = SPDFEPOCHUnixTime(unixtime, 'TOEPOCH', true);
%
%   % Convert a vector of unix times to CDF_EPOCH epochs.
%
%   epochs = SPDFEPOCHUnixTime(unixtimes, 'TOEPOCH', true);
%
%   Note: CDF_EPOCH does not have the time resolution for microseconds. The
%         conversion will not preserve it if the Unix time has that.

% HISTORY:
%   November 22, 2018  Mike Liu    The initial version.

%
% Process arguments.
%

if (nargin < 1)
    error('MATLAB:SPDFEPOCHUnixTime:inputArgumentCount', ...
          'SPDFEPOCHUnixTime requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFEPOCHUnixTime:outputArguments', ...
          'SPDFEPOCHUnixTime requires only one output argument.')
end

[args, msg] = parse_inputs(varargin{:});

if (~isempty(msg))
    error('MATLAB:spdfepochunixtime:badInputArguments', '%s', msg)
end

out = spdfepochunixtimec(time, args.toepoch);

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
    paramStrings = {'toepoch'};

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
         case 'toepoch'
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
