%
% Write (create/overwrite; open & close) an entire file from a 1D array of
% bytes. Cf irf.fs.read_file().
%
%
% ARGUMENTS
% =========
% filePath
%       Path to file to write.
% byteArray
%       uint8 column array.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-03.
%
function write_file(filePath, byteArray)
% PROPOSAL: Flag for permitting overwrite.
% NOTE: Code works if submitting string (char array), if it were not for the
%       assertion against it.
%   PROPOSAL: Permit char array.
%     CON: The caller can easily convert to uint8.
%     CON: There might be problems, subtleties relating to character tables,
%          Unicode, ASCII. Might not work as expected.

assert(isa(byteArray, 'uint8'), 'byteArray is not uint8.')
assert(iscolumn(byteArray), 'byteArray is not a column array.')



% fopen options:
% 'w'     open file for writing; discard existing contents
fileId = fopen(filePath, 'w');
% ~ASSERTION
if fileId == -1
  error('write_file:CanNotOpenFile', 'Can not open file: "%s"', filePath)
end

%===============
% Write to file
%===============
fwrite(fileId, byteArray);

errorCode = fclose(fileId);
% ~ASSERTION
if errorCode == -1
  error('write_file:CanNotCloseFile', 'Can not close file: "%s"', filePath)
end
end
