function epochObj = cdfepoch(varargin)
%CDFEPOCH Construct epoch object for CDF export.
%
%    E = CDFEPOCH(DATE) constructs a cdfepoch object where DATE is
%    a valid string (datestr) or number (datenum) representing a
%    date.  DATE may also be a cdfepoch object.
%
%    CDFEPOCH objects should be constructed to create EPOCH data in CDF's.
%    using SPDFCDFWRITE.  Note that a CDF epoch is the number of milliseconds
%    since 1-Jan-0000 and that MATLAB datenums are the number of days
%    since 0-Jan-0000.
%
%    See also SPDFCDFWRITE, DATENUM, SPDFCDFREAD, SPDFCDFINFO.

if (nargin == 0)
    s.date = [];
    epochObj = class(s, 'cdfepoch');
    return;
elseif (nargin > 1)
    error('MATLAB:cdfepoch:cdfepoch:tooManyInputs', ...
          'Too many input arguments.');
else
    input = varargin{1};
end

if isa(input,'cdfepoch')
    epochObj = input;
    return;
end

if iscellstr(input)
    input = char(input);
end

if ~ischar(input) & ~isnumeric(input)
    error('MATLAB:cdfepoch:cdfepoch:badInputs', ...
          'Input must be a number, string, cellstr, or cdfepoch object.');
end

% Initialize in case passed empty
s.date = [];

if ischar(input)
    % If the input is a string, then you have to convert
    cdfepochs = [];

    % Convert to MATLAB datenum.  If this bombs out, an invalid
    % datestr was passed to datenum.
    n = datenum(input);
else
    % It's numeric, so if it's a matrix, go element by element
    % and convert each and then reshape.
    n = input(:);
end

s = struct('date',num2cell((n - 1) * 24 * 3600000)');
s = s';

if isnumeric(input) & ~isempty(input)
    s = reshape(s, size(input));
end

epochObj = class(s, 'cdfepoch');

if isnumeric(input) & ~isempty(input)
    s = reshape(s, size(input));
end
