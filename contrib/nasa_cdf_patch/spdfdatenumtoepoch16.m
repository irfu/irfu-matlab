function epoch16 = spdfdatenumtoepoch16(varargin)
%SPDFDATENUMtoEPOCH16 Convert MATLAB's datenum to CDF_EPOCH16 values.
%
%    E = SPDFDATENUMtoEPOCH16(DATE) convert a DATE, a valid string (datestr) or
%        number (datenum) representing a date, to CDF_EPOCH16 value.
%
%    Note that a CDF_EPOCH16 is the number of picoseconds since 
%    1-Jan-0000 while MATLAB datenums are the number of days since 0-Jan-0000.
%
%    See also SPDFEPOCH16TODATENUM, SPDFCOMPUTEEPOCH16, SPDFENCODEEPOCH16.
%             SPDFPARSEEPOCH16, SPDFBREAKDOWNEPOCH16.

if (nargin > 1)
    error('MATLAB:SPDFDATENUMtoEPOCH16:SPDFDATENUMtoEPOCH16:tooManyInput', ...
          'Only one argument is allowed.');
elseif (nargin < 1)
    error('MATLAB:SPDFDATENUMtoEPOCH16:SPDFDATENUMtoEPOCH16:tooFewInput', ...
          'Only one argument is allowed.');
else
    input = varargin{1};
end

if iscellstr(input)
    input = char(input);
end

if ~ischar(input) & ~isnumeric(input)
    error('MATLAB:SPDFDATENUMtoEPOCH16:SPDFDATENUMtoEPOCH16:badInputs', ...
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

epoch16 = spdfdatenumtoepoch16c(n);

