function tt2000 = datenumtott2000(varargin)
%DATENUMtoTT2000 Convert MATLAB's datenum to CDF_TIME_TT2000 values.
%
%    E = DATENUMtoTT2000(DATE) convert a DATE, a valid string (datestr) or
%        number (datenum) representing a date, to CDF TT2000 value.
%
%    Note that a CDF epoch is the number of nanoseconds since 
%    1-Jan-0000T12:00:00 with leap seconds included and that MATLAB
%    datenums are the number of days since 0-Jan-0000.
%
%    See also CDFTT2000, TT2000TODATENUM, COMPUTETT2000, ENCODETT2000.
%             PARSETT2000.

%    binky
%    Copyright 2001-2006 The MathWorks, Inc.
%    $Revision: 1.1.1.1 $  $Date: 2012/05/17 14:33:09 $

if (nargin > 1)
    error('MATLAB:DATENUMtoTT2000:DATENUMtoTT2000:tooManyInput', ...
          'Only one argument is allowed.');
elseif (nargin < 1)
    error('MATLAB:DATENUMtoTT2000:DATENUMtoTT2000:tooFewInput', ...
          'Only one argument is allowed.');
else
    input = varargin{1};
end

if iscellstr(input)
    input = char(input);
end

if ~ischar(input) & ~isnumeric(input)
    error('MATLAB:DATENUMtoTT2000:DATENUMtoTT2000:badInputs', ...
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

tt2000 = datenumtott2000c(n);

