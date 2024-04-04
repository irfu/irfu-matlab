%
% Given a list of paths to files and directories, return a list of paths
% to files only. Directories are searched recursively for files. Also return
% dir() info for each path since it is a byproduct of the processing.
%
%
% ARGUMENTS
% =========
% fileDirPathsCa
%       Column cell array of paths to existing files and directories (mixed).
%
%
% RETURN VALUES
% =============
% filePathsCa
%       Column cell array of paths to files.
%       The (paths to) files in fileDirPathsCa are copied here.
%       The directories in fileDirPathsCa are searched recursively for files,.
% FsoiArray1
%       Column struct array.
%       dir() info for every path in filePathsCa.
%
%
% FSOI = File System Object Info
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [filePathsCa, FsoiArray1] = get_file_paths(fileDirPathsCa)
% NOTE: Function could be implemented recursively, but that is not needed since
%       dir() handles that.
%
% PROPOSAL: Redefine as generic, outside of BICAS.
%
% TODO-DEC: How handle non-existing files and directories?
%   NOTE: Only existing directories can be translated into files. Ambiguity
%         "only" exists for files.
%   PROPOSAL: Do not test for existing files.
%     PRO: Caller should test for that anyway.
%     CON: Can not return dir() info for files, even though the function does
%          produce it for directories anyway.
%   PROPOSAL: Test for (require) existing files.
%     PRO: Can return dir() info for all files.
%   PROPOSAL: Only return paths and dir() info for existing files. Also return
%             boolean array for input paths which are non-existent.
%
% PROPOSAL: Argument for specifying subset of files when recursing over
%           directories with dir().

assert(iscell(fileDirPathsCa) & iscolumn(fileDirPathsCa))



% Create empty struct array of the right size
% -------------------------------------------
% IMPLEMENTATION NOTE: Calling dir() just to make sure that the code uses
% the exact struct which it produces. This should make the code more
% future-proof.
FsoiEmptyArray = dir('~');
FsoiEmptyArray = FsoiEmptyArray([], 1);    % Column array.

bIsFile = arrayfun(@(pathCa) (exist(pathCa{1}, 'file') == 2), fileDirPathsCa, 'UniformOutput', true);
bIsDir  = arrayfun(@(pathCa) (exist(pathCa{1}, 'dir')  == 7), fileDirPathsCa, 'UniformOutput', true);

assert(all(bIsFile | bIsDir), 'Not all paths refer to existing files or directories.')

inputFilePathCa = fileDirPathsCa(bIsFile);
inputDirPathCa  = fileDirPathsCa(bIsDir);



% Pre-allocate
% ------------
% IMPLEMENTATION NOTE: A human caller may specify many explicit file paths using
% globbing in the OS. Therefore implemented function with pre-allocation of
% entries for file paths which in turn leads to splitting the processing of
% files and directories (though the splitting is not absolutely necessary).
fsoiCa = cell(1 + numel(inputFilePathCa), 1);
fsoiCa{1, 1} = FsoiEmptyArray;

for iFile = 1:numel(inputFilePathCa)
  fsoiCa{1+iFile, 1} = dir(inputFilePathCa{iFile});
end
FsoiArray1 = cat(1, fsoiCa{:});



for iDir = 1:numel(inputDirPathCa)
  path       = inputDirPathCa{iDir};
  FsoiArray2 = dir(fullfile(path, '**'));   % Recursive call to dir().
  FsoiArray2 = FsoiArray2(~[FsoiArray2.isdir]);
  FsoiArray2 = FsoiArray2(:);
  % CASE: FsoiArray2 is a column array.

  FsoiArray1 = [...
    FsoiArray1;
    FsoiArray2];
end

filePathsCa = arrayfun(@(Oi) (fullfile(Oi.folder, Oi.name)), FsoiArray1, 'UniformOutput', false);



assert(iscolumn(filePathsCa))
assert(iscolumn(FsoiArray1))
assert(numel(filePathsCa) == numel(FsoiArray1))

end
