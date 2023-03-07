function tt2000Obj = cdftt2000(varargin)
%CDFTT2000 Construct a cdftt2000 object for CDF export.
%
%    E = CDFTT2000(DATE) constructs a cdftt2000 object where DATE is
%    a valid UTC string (datestr) or number (mxINT64_CLASS) representing a
%    TT2000 date.  DATE may also be a cdftt2000 object.
%
%    CDFTT2000 objects should be constructed to create TT2000 data in CDF's
%    using CDFWRITE.  Note that a CDF TT2000 is the number of nanoseconds
%    since 1-Jan-2000 at 12:00:00.
%
%    See also SPDFCDFWRITE, SPDFCDFREAD, SPDFCDFINFO, SPDFENCODETT2000,
%             SPDFCOMPUTETT2000, SPDFPARSETT2000, SPDFBREAKDOWNTT2000,
%             SPDFCDFLEAPSECONDSINFO.

%    $Revision: 1.1.1.1 $  $Date: 2023/01/24 16:42:29 $

if (nargin == 0)
    s.date = [];
    tt2000Obj = class(s, 'cdftt2000');
    return;
elseif (nargin > 1)
    error('MATLAB:cdftt2000:cdftt2000:tooManyInputs', ...
          'Too many input arguments.');
else
    input = varargin{1};
end

if isa(input,'cdftt2000')
    tt2000Obj = input;
    return;
end

if iscellstr(input)
    input = char(input);
end

if ~ischar(input) & ~isnumeric(input) & ~iscell(input)
    error('MATLAB:cdftt2000:cdftt2000:badInputs', ...
          'Input must be a number, string, cellstr, or cdftt2000 object.');
end

% Initialize in case passed empty
s.date = [];

if ischar(input)

    % Convert to TT2000 values. 
    n = spdfparsett2000(input);
else
    % It's numeric, so if it's a matrix, go element by element
    % and convert each and then reshape.
    if iscell(input)
      if (isa(input, 'int64'))
        n = int64(input{:});
      else
        n = input{:};
      end
    else
      if (isa(input, 'int64'))
        n = int64(input(:));
      else
        n = input(:);
      end
    end
end

s = struct('date',n);
s = s';
if isnumeric(input) & ~isempty(input)
    s = reshape(s, size(input));
end

tt2000Obj = class(s, 'cdftt2000');
if isnumeric(input) & ~isempty(input)
    s = reshape(s, size(input));
end

end
