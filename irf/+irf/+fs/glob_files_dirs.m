%
% Select files and directories in a globbing-like way. Uses regex patterns
% instead of globbing syntax.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% rootDirPath
%       Directory path. Selected files/directories will be located under this
%       directory (recursively).
% regexPatternList
%       1D cell array of regex patterns (strings). One component for every part
%       of the relative paths (i.e. excluding the root directory name itself). A
%       regex matching ANY SUBSTRING counts as a match. To match only the entire
%       object name, use regex commands ^ and $.
% varargin
%       Settings struct and/or list of pairs of arguments representing settings
%       (string key + value) to override the defaults. Values are passed on to
%       irf.fs.recurse_directory_tree(). See implementation for
%       meaning and default values.
%
%
% RETURN VALUE
% ============
% objectInfoList
%       Struct array. Can be zero-length.
%         .name
%         .date
%         .bytes
%         .isdir
%         .datenum
%         .relativePath
%         .recursionDepth
%         .fullPath
%       NOTE: If one submits an empty cell array of regex patterns, then it is a
%       one-component struct array for the root directory (i.e. it generalizes
%       nicely).
%       NOTE: Some fields have field names inherited from the dir() command.
%
%
% NOTES
% =====
% MATLAB does not appear to have any similar functionality. "glob" in the
% function name is misleading since it really uses regular expressions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-04-06
%
function objectInfoList = glob_files_dirs(rootDirPath, regexPatternList, varargin)
%
% PROPOSAL: Argument for choosing types of objects: files, directories, files & directories.
%   PROPOSAL: Implement as setting
%       NOTE: Settings are passed on to irf.fs.recurse_directory_tree. Problem?
%
% PROPOSAL: Change function name to not use "glob".
%   PRO: "glob" is deceiving when uses regex.
%       CON: Wikipedia hints that it is OK, that "globbing" refers to matching files (and paths):
%            "In computer programming, glob patterns specify sets of filenames with wildcard characters."
%   PROPOSAL: select_files, specify_files, *_by_name, *_by_path
%   PROPOSAL: regexp_files_dirs
%   PROPOSAL: Something with "path".
% PROPOSAL: Not follow naming convention from dir().
%
% PROPOSAL: Option for specifying regexp for entire (relative) path, including slash.
%   TODO-NI: Path details become important. Paths ending/not ending with slash; current directory, parent directory.
%   PROPOSAL: If regexPatternList is a string (not cell array of strings).


% ASSERTION
assert(iscell(regexPatternList), 'regexPatternList is not a cell array.')

objectInfoList = irf.fs.recurse_directory_tree(...
  rootDirPath, ...
  @(args) FileFunc         (args, regexPatternList), ...
  @(args) DirFunc          (args, regexPatternList), ...
  @(args) ShouldRecurseFunc(args, regexPatternList), ...
  varargin{:});

if isempty(objectInfoList)
  objectInfoList = irf.ds.empty_struct([0,1], ...
    'name', 'date', 'bytes', 'isdir', 'datenum', ...
    'relativePath', 'recursionDepth');
end

end



function result = FileFunc(args, regexPatternList)

if (args.recursionDepth == length(regexPatternList)) ...
    && matchesRegex(args.dirCmdResult.name, regexPatternList{end})

  result = args.dirCmdResult;
  result.recursionDepth = args.recursionDepth;
  result.relativePath   = args.relativePath;
  result.fullPath       = args.fullPath;
else
  result = [];
end

end



function result = DirFunc(args, regexPatternList)

result = [];

%fprintf('args.relativePath = %s   ; args.hasRecursedOverChildren = %i\n', args.relativePath, args.hasRecursedOverChildren);
if args.hasRecursedOverChildren ...
    && (args.recursionDepth == length(regexPatternList))
  % Add current directory to the directory results.
  result = args.dirCmdResult;
  result.recursionDepth = args.recursionDepth;
  result.relativePath   = args.relativePath;
  result.fullPath       = args.fullPath;
end

% Add childrens' results to the directory results.
for i = 1:length(args.childrenResultsList)
  result = [result; args.childrenResultsList{i}];
end

end



function shouldRecurse = ShouldRecurseFunc(args, regexPatternList)

if args.recursionDepth == 0
  shouldRecurse = true;
elseif args.recursionDepth > length(regexPatternList)
  shouldRecurse = false;
else
  regexPattern = regexPatternList{args.recursionDepth};
  shouldRecurse = matchesRegex(args.dirCmdResult.name, regexPattern);
end
%fprintf('args.relativePath = %s   ; args.shouldRecurse = %i\n', args.relativePath, shouldRecurse);

end



function matches = matchesRegex(str, regexPattern)
matches = ~isempty(regexp(str, regexPattern));
end
