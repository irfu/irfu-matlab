function out = SPDFENCODEEPOCH(epoch, varargin)
%SPDFENCODEEPOCH encodes an epoch of CDF_EPOCH data type, a double value or
%cdfepoch object..
%
%   OUT = SPDFENCODEEPOCH(epoch) 
%         Returns a UTC string.
%
%     epoch                An epoch
%
%   The encoded epoch string will have the following ISO 8601 format:
%       yyyy-mm-ddThh:mm:ss.mmm, e.g., "2000-01-01T12:34:56.123"
%   Originally, it was in this form:
%       dd-mmm-yyyy hh:mm:ss.mmm, e.g., "01-Jan-2000 12:34:56.123"
%
%   OUT = SPDFENCODEEPOCH(epoch, 'Format', FORMAT) encodes the UTC
%   string into the specified format. FORMAT is a number from 0 to 4.
%   FORMAT:
%     0: dd-mmm-yyyy hh:mm:ss.mmm, e.g., "01-JAN-2000 12:34:56.123"
%     1: yyyymmdd.ddddddd, e.g., "20000101.1200000"
%     2: yyyymmddhhmmss, e.g., "20000101123456"
%     3: yyyy-mm-ddThh:mm:ss.mmmZ, e.g., "2000-01-01T12:34:56.123Z"
%     4: yyyy-mm-ddThh:mm:ss.mmm, e.g., "2000-01-01T12:34:56.123"
%   where mmm is milliseconds.
%   Format: 0 is the only allowed form for cdfepoch object.
%
%   Note: If the epoch values come from spdfcdfread function call, the values
%         can be in one of the three forms: in cdfepoch object, in MATLAB
%         datenum (the default), or in their original CDF_EPOCH 
%         data via an extra 'keepepochasis' option. This function works for
%         cdfepoch objects or the data retrieved with 'keepepochasis' option.
%         For datenum values, use MATLAB's datestr instead.
%
%   Examples:
%
%   % Encode epoch from date/time: 2012-10-10T10:10:10.010Z:
%
%   utc = '2012-10-10T10:10:10.010';
%   epoch = UTC2CDFEpoch(utc);
%   SPDFENCODEEPOCH(epoch)
%   ans =
%       '10-Oct-2012 10:10:10.010'
%   SPDFENCODEEPOCH(epoch, 'format', 3) 
%   ans =
%       '2012-10-10T10:10:10.010Z'
%
%   % Acquire 'Epoch' variable data as is (double values) from 'sample' CDF
%   % and encode the epoch values.
%
%   epochs = spdfcdfread('Sample', 'Variables', {'Epoch'}, 'KEEPEPOCHASIS', true);
%   spdfencodeepoch(epochs)
%
%   See also CDFEPOCH, SPDFBREAKDOWNEPOCH, SPDFCOMPUTEEPOCH, SPDFPARSEEPOCH

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.
%   October 16, 2018  Mike Liu   The default format is now 4 (from 0).

if (nargin < 1)
    error('MATLAB:SPDFENCODEEPOCH:inputArgumentCount', ...
          'SPDFENCODEEPOCH requires at least one input argument.')
end

if (nargout > 1)
    error('MATLAB:SPDFENCODEEPOCH:outputArguments', ...
          'SPDFENCODEEPOCH requires only one output argument.')
end

[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:SPDFENCODEEPOCH:badInputArguments', '%s', msg)
end
if (isa(epoch,'cdfepoch'))
  if (length(epoch) > 1)
    for p = 1:length(epoch)
      if (~isa(epoch, 'cell'))
        dataaa(p,1) = datestr((todatenum(epoch(p,1))), 0);
      else
        if (length(epoch{p}) > 1)
          for q = 1:length(epoch{p})
            dataaa(q) = datastr((todatenum(epoch{p}(q,1))), 0);
          end
        else
          dataaa(p,1) = datastr((todatenum(epoch{p})), 0);
        end
      end
    end
    out = dataaa;
  end
else
  out = spdfencodeepochc(epoch, args.Format);
end
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

