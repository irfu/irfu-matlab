%
% Write (create/overwrite; open & close) an entire file from a 1D array of
% bytes. Cf read_file.
%
%
% ARGUMENTS
% =========
% byteArray
%       uint8 1D array.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-03.
%
function write_file(filePath, byteArray)
% PROPOSAL: Flag for permitting overwrite.

assert(isa(byteArray, 'uint8'), 'byteArray is not uint8.')
irf.assert.vector(byteArray)

% fopen options:
% 'w'     open file for writing; discard existing contents
fileId = fopen(filePath, 'w');
% ~ASSERTION
if fileId == -1
  error('write_file:CanNotOpenFile', 'Can not open file: "%s"', filePath)
end

fwrite(fileId, byteArray);


errorCode = fclose(fileId);
% ~ASSERTION
if errorCode == -1
  error('write_file:CanNotCloseFile', 'Can not close file: "%s"', filePath)
end
end
