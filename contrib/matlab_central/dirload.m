function out = dirload(file,varargin)
%DIRLOAD Like load, but each var in its own file
if nargin < 1 || isempty(file); file = 'matlab'; end
varNames = varargin;
if ~exist(file, 'dir')
    error('Not a dirsave dir: %s', file);
end
if isempty(varNames)
    d = dir(file);
    varNames = regexprep(setdiff(ls(file), {'.','..'}), '\.mat$', '');
end

out = struct;
for i = 1:numel(varNames)
    name = varNames{i};
    tmp = load(fullfile(file, [name '.mat']));
    out.(name) = tmp.(name);
end

if nargout == 0
    for i = 1:numel(varNames)
        assignin('caller', varNames{i}, out.(varNames{i}));
    end
    clear out
end
