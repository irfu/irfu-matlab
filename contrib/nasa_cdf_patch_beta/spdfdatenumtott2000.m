function tt2000 = spdfdatenumtott2000(varargin)
%SPDFDATENUMtoTT2000 Convert MATLAB's datenum to CDF_TIME_TT2000 values.
%
%    E = SPDFDATENUMtoTT2000(DATE) convert a DATE, a valid string (datestr) or
%        number (datenum) representing a date, to CDF TT2001 value.
%
%    Note that a CDF TT2000 is the number of nanoseconds since 
%    1-Jan-0000T12:00:00 with leap seconds included and that MATLAB
%    datenums are the number of days since 0-Jan-0000.
%
%    See also CDFTT2000, SPDFTT2000TODATENUM, SPDFCOMPUTETT2000, SPDFENCODETT2000,
%             SPDFPARSETT2000, SPDFBREAKDOWNTT2000.


if (nargin > 1)
    error('MATLAB:SPDFDATENUMtoTT2000:SPDFDATENUMtoTT2000:tooManyInput', ...
          'Only one argument is allowed.');
elseif (nargin < 1)
    error('MATLAB:SPDFDATENUMtoTT2000:SPDFDATENUMtoTT2000:tooFewInput', ...
          'Only one argument is allowed.');
else
    input = varargin{1};
end

if iscellstr(input)
    input = char(input);
end

if ~ischar(input) & ~isnumeric(input)
    error('MATLAB:SPDFDATENUMtoTT2000:SPDFDATENUMtoTT2000:badInputs', ...
          'Input must be a number, string, cellstr.');
end

if ischar(input)
    % If the input is a string, then you have to convert

    % Convert to MATLAB datenum.  If this bombs out, an invalid
    % datestr was passed to datenum.
    n = datenum(input);
else
    % It's numeric, so if it's a matrix, go element by element
    % and convert each and then reshape.
    n = input(:);
end

tt2000 = spdfdatenumtott2000c(n);

