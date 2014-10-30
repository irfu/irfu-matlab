function SPDFCDFLEAPSECONDSINFO(varargin)
%SPDFCDFLEAPSECONDSINFO shows the information how the leap seconds is used by
%CDF.
%
%   SPDFCDFLEAPSECONDSINFO() displays the basic setup the CDF library uses to 
%   acquire the leap second table.
%
%   SPDFCDFLEAPSECONDSINFO('DUMP', TF) the 'DUMP' option will show the contents
%   of the leap second table, in addition to the basic information, if TF is
%   true.
%
%   See also CDFTT2000, SPDFENCODETT2000, SPDFCOMPUTETT2000, SPDFPARSETT2000,
%            SPDFBREAKDOWNTT2000.

% HISTORY:
%   August 16, 2011  Mike Liu    The initial version.

%
% Process arguments.
%

[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:SPDFCDFLEAPSECONDSINFO:badInputArguments', '%s', msg)
end

spdfcdfleapsecondsinfoc(args.Dump);

%%%
%%% Function parse_inputs
%%%

function [args, msg] = parse_inputs(varargin)
% Set default values
args.Dump = false;
msg = '';
% Parse arguments based on their number.
if (nargin > 0)
    paramStrings = {'dump'};

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
       case 'dump'

           if (k == length(varargin))
               msg = 'No dump specified.';
               return
           else

               dump = varargin{k + 1};
               if (numel(dump) ~= 1)
                   msg = 'The "Dump" value must be a scalar logical.';
               end

               if (islogical(dump))
                   args.Dump = dump;
               elseif (isnumeric(dump))
                   args.Dump = logical(dump);
               else
                   msg = 'The "Dump" value must be a scalar logical.';
               end
           end

       end  % switch
    end  % for 
                   
end  % if (nargin > 1)


