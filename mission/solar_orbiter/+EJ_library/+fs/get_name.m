%
% Extract file/directory/object name from path.
% 
% NOTE: Will interpret anything after last filesep as filename, including ".." or ".", or empty string for paths that
% end with filesep.
% NOTE: Pure string operation. Path is not required to exist.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created circa 2020-04-29.
%
function name = get_name(path)
    % PROPOSAL: Generalize to accepting cell array of paths.
    %   CON: Easy to use cellfun.
    %       Ex: cellfun(@EJ_library.fs.get_name, {DsmdArray1.path}', 'UniformOutput', false)
    %   CON: Provides no extra value, e.g. speed, except for convenience over using cellfun.
    %
    % PROPOSAL: Return parentPath too (second return argument for compatibility).
    %   PROPOSAL: Change name get_name_parent.
    
    [~, basename, suffix] = fileparts(path);
    name = [basename, suffix];
end
