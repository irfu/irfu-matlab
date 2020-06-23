%
% Generic function for recursing over the files & directories under an arbitrary directory.
% 
% Recurses over all, or a subset of, files and directories under a directory subtree ("the root path/directory"). Calls
% one arbitrary function for every directory and another arbitrary function for every file. A third arbitrary function
% determines whether the children of any directory should be recursed over (useful for improving speed by avoiding
% processing).
%
%
%
% ARGUMENTS
% =========
% rootDirPath
% FileFunc          : Function pointer that is called for every file that is recursed over.
%                           result = file_func(argStruct).
%                     result    : An arbitrary value that is passed on to the call for "DirFunc" for the parent
%                                 directory, via the "childrenResultsList" argument.
%                     argStruct : Struct with fields
%                          .rootDirPath
%                          .relativePath
%                          .fullPath
%                          .dirCmdResult
%                          .recursionDepth
%
% DirFunc           : Function pointer, that is called for every directory recursed over, including the root directory
%                     ("rootDirPath"). It is implicit that this function should not need to recurse over the directory
%                     tree itself.
%                           result = dir_func(ArgStruct)
%                     result    : An arbitrary value passed on to the call for "DirFunc" for the parent directory. For
%                                 the root directory, the value is returned by "recurse_directory_tree" itself.
%                     ArgStruct : Struct with fields
%                          .rootDirPath
%                          .relativePath
%                          .fullPath
%                          .dirCmdResult
%                          .recursionDepth
%                          .childrenResultsList     : Cell array of "result" values returned from FileFunc and DirFunc
%                                                     when called for the directory's children.
%                          .hasRecursedOverChildren : The result of ShouldRecurseFunc for the same directory.
%                                                     Is used to distinguish the absence of children from the absence of
%                                                     recursion over children.
%
% ShouldRecurseFunc : Function pointer:
%                           shouldRecurse = should_recurse(ArgStruct)
%                     shouldRecurse : True/false. True iff "recurse_directory_tree" should recurse over the children of
%                                     the specified directory.
%                                     NOTE: DirFunc will always be called for the directory specified here.
%                     ArgStruct     : Struct with fields
%                          .rootDirPath
%                          .relativePath
%                          .fullPath
%                          .dirCmdResult
%                          .recursionDepth
%                     
% varargin          : Optional settings on format determined by EJ_library.utils.interpret_settings_args .
%                     useRelativeDirectorySlash : True/false. Whether "relativePath" for directories should have a trailing
%                                                 slash/backslash or not.
%                     useRootRelativePathPeriod : True/false. Whether ".relativePath" for rootDirPath (the "root path")
%                                                 itself should be represented by a period (instead of an empty string).
%
%
%
% RETURN VALUE
% ============
% result            : What DirFunc returns for the root path. Its format and content are determined by FileFunc, DirFunc
%                     and ShouldRecurseFunc.
%
%
%
% VARIABLE NAMING CONVENTION
% ==========================
% rootDirPath    : Path to "root path/directory", i.e. the root of the subtree to be recursed over.
% relativePath   : Path to a file/directory relative to "rootDirPath".
%                  NOTE: Code uses an empty string instead of "." (or "./") to represent the root path.
% fullPath       : Path to a file/directory beginning with "rootDirPath". Is thus still relative if "rootDirPath" is.
% dirCmdResult   : The struct with information for this object that is returned by MATLAB's "dir", 
%                  i.e. NOT information on all the children, only the object itself.
% recursionDepth : Number of levels into the rootDirPath subdirectory tree. rootDirPath is level zero.
%                  Files and directories have the same level.
%
%
%
% NOTE: Function might be used with Lapdog. Should therefore be compatible with MATLAB R2009a.   % How? Why?
%
% 
%
% DESIGN RATIONALE / IMPLEMENTATION NOTE
% ======================================
% The purpose of this function is partly to make it easy to implement many other generic functions.
% Examples:
% (1) Constructing a list of files and/or directories under a directory path, incl. selectively.
% (2) Summing up the file sizes of files in a directory tree.
% (3) Constructing a recursive data structure analogous to, or based upon, the subtree.
% (4) Executing code for every file and/or directory in a subtree.
%     Alt 1) If needes shared state, collect list of paths to file and/or directories. Then iterate over them.
%     Alt 2) Execute code in FileFunc and/or DirFunc.
%
% IMPLEMENTATION NOTE: Function pointers take structs as arguments instead of conventional lists of arguments to,
% (1) make it easier(?) to write one-line anonymous functions,
% (2) make it easier to add new (effective) arguments to the function pointers while maintaining backward compatibility
%     with code that uses recurse_directory_tree,
% (3) avoid having to remember the order of the arguments when specifying function pointers (both for internal
%     implementation and for external callers),
% (4) make it possible to use the same function pointer for both FileFunc and DirFunc (can use "dirCmdResult.isDir" to
%     distinguish calls to FileFunc and DirFunc). Therefore, argStruct has as many identical field names as possible.
%
%
%
% HANDLING OF SPECIAL CASES -- NEEDS TO BE CONFIRMED
% ==================================================
% Should work for:
%    (1) rootDirPath = "/" (Linux).
%    (2) rootDirPath = Path without the name for the final subdirectory, e.g. ".", "/qwe/asd/.", "/qwe/asd/.."
%    (3) Empty directory (or subdirectory)
% Not tested for: 
%    (1) Symbolic links, hard links, devices etc but it should be possible to derive the behaviour for symbolic links 
%        using the notes below.
%
%
%
% NOTE: "dir" COMMAND RESULTS (EXAMPLES) FOR SYMBOLIC LINKS
% =========================================================
% Symbolic link to existing file:
%       name: 'filelink'
%       date: '28-feb-2017 17:38:21'
%      bytes: 5
%      isdir: 0
%    datenum: 7.3675e+05
% Symbolic link to existing directory:
%       name: 'dirlink'
%       date: '28-feb-2017 17:38:21'
%      bytes: 0
%      isdir: 1
%    datenum: 7.3675e+05
% Symbolic link to nothing:
%       name: 'badlink'
%       date: ''
%      bytes: []
%      isdir: 0
%    datenum: []
%
%
%
% Initially created by Erik P G Johansson, 2016-11-29.
%
function result = recurse_directory_tree(rootDirPath, FileFunc, DirFunc, ShouldRecurseFunc, varargin)
% GUIDING NEEDS / RATIONALE / MOTIVATION
% ======================================
% NEED: Få namn för varje objekt.
% NEED: Relativa sökvägar för varje objekt.
% NEED: Få tillgång till (1D) lista med objekt i katalogträd.
% NEED: "Snabb" kod för stora träd, iaf om bara behöver mindre delar av det.
%     NEED: Utesluta delar av katalogträd för att snabba upp iteration.
%         NEED: Utesluta rekursion under vissa kataloger.
%         NEED: Utesluta rekursion bort visst rekursionsdjup.
% NEED: Ska kunna användas för att lätt implementera
%     (1) Ta fram objektlista mha globbing (fixt rekursionsdjup; ~regex).
%     (2) Producera lista med alla filer, eller filtrerad delmängd därav, under katalog.
%     (3) Sökning efter fil(er) (baserat på filnamn/plats, innehåll).
% NEED: Lätt kunna ta fram aggregerande information rekursivt för helt träd.
%     Ex: Total storlek (bytes).
%     Ex: Antal filer/kataloger.
%     Ex: Aggregerade data från flera filer.
%     (Ex: Lista med objekt (eller urval därav).)
%     ==> Exekvera kod för varje objekt: fil-->data, katalog+barndata-->data

% PROPOSAL: Argument för maximalt rekursionsdjup.
%   CON: Behövs inte. ShouldRecurseFunc kan lätt implementera det om det behövs.
% PROPOSAL: Tillåt kombination  Settings.useRelativeDirectorySlash && ~Settings.useRootRelativePathPeriod.
%           Lägg då INTE till snedstreck för rootkatalogen.
% PROPOSAL: Avskaffa useRelativeDirectorySlash. Aldrig sluta med slash.
%   NOTE: Se MATLAB-anteckningsfil.
% PROPOSAL: Merge ShouldRecurseFunc into DirFunc.
%   ~CON: Not obvious how this works with algorithm.

    

    % ASSERTION
    if ~exist(rootDirPath, 'dir')
        error('Can not find directory "%s".', rootDirPath)
    end


    
    % Set default settings.
    DEFAULT_SETTINGS.useRelativeDirectorySlash = false;
    DEFAULT_SETTINGS.useRootRelativePathPeriod = false;
    Settings = EJ_library.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);

    % ASSERTION
    if Settings.useRelativeDirectorySlash && ~Settings.useRootRelativePathPeriod
        error('Illegal combination of settings.')
    end
    
    % Convert settings into values to actually use (set once here for speed; to avoid repetition).
    if Settings.useRelativeDirectorySlash ; relativeDirectoryPathSuffix = '/';
    else                                    relativeDirectoryPathSuffix = '';
    end
    if Settings.useRootRelativePathPeriod ; relativePathRoot = '.';    % Should NOT include trailing slash. Is added automatically if needed.
    else                                    relativePathRoot = '';     % Should NOT include trailing slash. Is added automatically if needed.
    end



    dirCmdResultRootPath = get_dir_cmd_result_for_single_object(rootDirPath);

    % NOTE: Using empty string to represent relative path to rootDirPath instead of period. The string is used for is
    % used for building other relative paths (with "fullfile"). If one uses '.', or './' then one gets "ugly" relative
    % paths beginning with "./" which is unnecessary for the actual children of the directory.
    result = recurse_directory_tree_INTERNAL(...
        rootDirPath, '', ...
        dirCmdResultRootPath, 0, ...
        FileFunc, DirFunc, ShouldRecurseFunc, ...
        relativeDirectoryPathSuffix, relativePathRoot);
end



% Function that should be called for every subdirectory, and never for files.
%
% ARGUMENTS
% =========
% rootDirPath         : The path to the root of the ENTIRE recursion. This value stays constant for all calls.
% relativePath        : Path to directory relative to "rootDirPath".
%                       NOTE: Uses empty string to represent "rootDirPath".
%                       NOTE: Never has trailing slash/backslash.
% dirCmdResultCurrent : The "dir" command result (data, record) that pertains to this specific directory, i.e. NOT
%                       the data/records for all of its immediate children (files and directories) contained in it.
% recursionDepth      : Zero for the first call, i.e. for rootDirPath (relativePath empty). Incremented by one for every
%                       step into the directory structure.
% FileFunc, DirFunc, ShouldRecurseFunc : Same as arguments to main function.
% 
% RETURN VALUE
% ============
% result              : [] (not cell) if ShouldRecurseFunc returns false for this directory.
% 
function result = recurse_directory_tree_INTERNAL(...
        rootDirPath, relativePath, ...
        dirCmdResultCurrent, recursionDepth, ...
        FileFunc, DirFunc, ShouldRecurseFunc, ...
        relativeDirectoryPathSuffix, relativePathRoot)
    
    % "dir" command results to remove/ignore, as identified by their .name fields.
    NON_CHILDREN_NAMES = {'.', '..'};
    
    % Value to return when shouldRecurse == false
    NO_RECURSE_DIRECTORY_RESULT = {};

    %fprintf(1, 'BEGIN: recurse_directory_tree_INTERNAL("%s")\n', relativePath)
    
    
    % Set "relativeDirPath". Only used for calling ShouldRecurseFunc and DirFunc.
    if recursionDepth == 0
        relativeDirPath = [relativePathRoot, relativeDirectoryPathSuffix];
    else
        relativeDirPath = [relativePath,     relativeDirectoryPathSuffix];
    end
    fullDirPath = fullfile(rootDirPath, relativeDirPath);
    
    
    
    %========================
    % Call ShouldRecurseFunc
    %========================
    willRecurseOverChildren = ShouldRecurseFunc(struct(...
        'rootDirPath',    rootDirPath, ...
        'relativePath',   relativeDirPath, ...
        'fullPath',       fullDirPath, ...
        'dirCmdResult',   dirCmdResultCurrent, ...
        'recursionDepth', recursionDepth));

    %=============================================
    % Derive "childrenResultsList" is supposed to.
    %=============================================
    if willRecurseOverChildren

        %==========================================================
        % Get information on the children of the current directory
        %==========================================================
        currentFullPath   = fullfile(rootDirPath, relativePath);
        dirCmdResultsList = dir(currentFullPath);                 % Returns column vector of structs.
        % Remove non-children from dir command results.
        iDelete = ismember({dirCmdResultsList.name}, NON_CHILDREN_NAMES);
        dirCmdResultsList(iDelete) = [];

        %=======================================
        % Iterate over the directory's children
        %=======================================
        childrenResultsList  = cell(length(dirCmdResultsList), 1);
        for iChild = 1:length(dirCmdResultsList)
            dirCmdResultChild = dirCmdResultsList(iChild);
            childRelativePath = fullfile(relativePath, dirCmdResultChild.name);

            if dirCmdResultChild.isdir
                % CASE: Child is a directory.

                %================
                % RECURSIVE CALL
                %================
                childrenResultsList{iChild} = recurse_directory_tree_INTERNAL(...
                    rootDirPath, childRelativePath, ...
                    dirCmdResultChild, recursionDepth+1, ...
                    FileFunc, DirFunc, ShouldRecurseFunc, ...
                    relativeDirectoryPathSuffix, relativePathRoot);
            else
                % CASE: Child is a file.

                %===============
                % Call FileFunc
                %===============
                %fprintf(1, 'Calling FileFunc("%s")\n', childRelativePath)
                %childrenResultsList{iChild} = FileFunc(rootDirPath, childRelativePath, dirCmdResultChild);
                childrenResultsList{iChild} = FileFunc(struct(...
                    'rootDirPath',    rootDirPath, ...
                    'relativePath',   childRelativePath, ...
                    'fullPath',       fullfile(rootDirPath, childRelativePath), ...
                    'dirCmdResult',   dirCmdResultChild, ...
                    'recursionDepth', recursionDepth+1));
            end
        end

    else
        childrenResultsList = NO_RECURSE_DIRECTORY_RESULT;
    end

    %fprintf(1, 'Calling DirFunc("%s")\n', relativePath)
    
    %==============
    % Call DirFunc
    %==============
    % NOTE: Requires cell braces around "childrenResultsList" to prevent "struct" from creating an analogous array of
    % structs. (Other fields are identical for all array components.)
    result = DirFunc(struct(...
        'rootDirPath',              rootDirPath, ...
        'relativePath',             relativeDirPath, ...
        'fullPath',                 fullDirPath, ...
        'dirCmdResult',             dirCmdResultCurrent, ...
        'recursionDepth',           recursionDepth, ...
        'childrenResultsList',      {childrenResultsList}, ...
        'hasRecursedOverChildren',  willRecurseOverChildren));

    %fprintf(1, 'END:   recurse_directory_tree_INTERNAL("%s")\n', relativePath)
end



% Return information for specific directory using "dir" (i.e. NOT information for all the children).
% Utility function for clarifying the code.
%
% Should be called exactly once, for the root directory.
function dirCmdResult = get_dir_cmd_result_for_single_object(path)
    % TODO-DECISION: How handle "/"?
    
    % IMPLEMENTATION NOTE: Need correct path to later determine the name of the directory (e.g. for relative paths, "/", "..").
    % EJ_library.fs.get_abs_path might not be enough.
    absPath = EJ_library.fs.get_abs_path(path);   
    
    dirCmdResultsList = dir(absPath);

    % Replace "." (current directory) with the actual name of the current directory.
    [junk, baseName, ext] = fileparts(absPath);   % IMPLEMENTATION NOTE: Keep compatible with MATLAB R2009a ==> Do NOT use "~" notation.
    iPath = find(strcmp({dirCmdResultsList.name}, '.'));
    assert(isscalar(iPath), 'Can not find exactly one instance of object named "." .')
    dirCmdResult = dirCmdResultsList(iPath);
    
    % IMPLEMENTATION NOTE: Empirically, "dir" returns zero entries for non-readable directories WITHOUT THROWING ANY
    % EXCEPTION. This can give very non-intuitive errors.
    assert(length(dirCmdResult) == 1, '"dir" returned zero entries for "%s". Non-existent? No read permissions?', absPath)
    dirCmdResult.name = [baseName, ext];            
end

