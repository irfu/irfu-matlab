%
% Read an entire file as a 1D array of byte values.
%
%
% RETURN VALUE
% ============
% byteArray
%       Column array of uint8.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-10.
%
function byteArray = read_file(filePath)
fileId = fopen(filePath);

% ASSERTION
if fileId == -1
  error('read_file:CanNotOpenFile', 'Can not open file: "%s"', filePath)
end

doubleArray = fread(fileId);
fclose(fileId);

% IMPLEMENTATION NOTE: Can not use int8() function directly since it returns
% signed int8 and the function saturates, i.e. e.g. int8(200)==int8(127).

% ASSERTIONS: Assert that there is no misunderstanding of what is returned from fread.
assert(all(doubleArray >= 0))
assert(all(doubleArray <=255))

byteArray = uint8(doubleArray);
assert(all(doubleArray == byteArray))

% Normalize to column array. (Only needed for empty array).
byteArray = byteArray(:);
end
