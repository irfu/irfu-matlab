function epoch = spdfdatenumtoepoch(varargin)
%SPDFDATENUMtoEPOCH Convert MATLAB's datenum to CDF_EPOCH values.
%
%    E = SPDFDATENUMtoEPOCH(DATE) convert a DATE, a valid string (datestr) or
%        number (datenum) representing a date, to CDF_EPOCH value.
%
%    Note that a CDF epoch is the number of milliseconds since 
%    1-Jan-0000 while MATLAB datenums are the number of days since 0-Jan-0000.
%
%    See also CDFEPOCH, SPDFEPOCHTODATENUM, SPDFCOMPUTEEPOCH, SPDFENCODEEPOCH.
%             SPDFPARSEEPOCH, SPDFBREAKDOWNEPOCH.

if (nargin > 1)
    error('MATLAB:SPDFDATENUMtoEPOCH:SPDFDATENUMtoEPOCH:tooManyInput', ...
          'Only one argument is allowed.');
elseif (nargin < 1)
    error('MATLAB:SPDFDATENUMtoEPOCH:SPDFDATENUMtoEPOCH:tooFewInput', ...
          'Only one argument is allowed.');
else
    input = varargin{1};
end

if iscellstr(input)
    input = char(input);
end

if ~ischar(input) & ~isnumeric(input)
    error('MATLAB:SPDFDATENUMtoEPOCH:SPDFDATENUMtoEPOCH:badInputs', ...
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

epoch = spdfdatenumtoepochc(n);

