function dirdelete(file,varargin)
%DIRDELETE Delete variables from a dirsave dir
if nargin < 1 || isempty(file); file = 'matlab'; end
varNames = varargin;
for i = 1:numel(varNames)
    delete(fullfile(file, [varNames{i} '.mat']));
end
