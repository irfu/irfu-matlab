function out = dirwhos(file,varargin)
%DIRWHOS List variable names in a dirsave dir
%
%  DIRWHOS(dirname) list variables in dirsave directory dirname
%
%  DIRWHOS(dirname,varname) list variable varname
%  OUT = DIRWHOS(dirname,varname) return output, same as whos

% developed from http://stackoverflow.com/questions/4268044/deleting-variables-from-a-mat-file
% SPDX-License-Identifier: CC-BY-SA-2.5

if nargin < 1 || isempty(file); file = 'matlab'; end
if nargin == 1
  out = regexprep(setdiff(ls(file), {'.','..'}), '\.mat$', '');
  out=out{1};
end
if nargin == 2
  fileName = varargin{1};
  filePath =[file filesep fileName '.mat'];
  if exist(filePath,'file')
    out = whos('-file',filePath,fileName);
  else
    out = [];
  end
end
