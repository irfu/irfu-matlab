%
% Wrapper around bicas.tools.batch.main() intended to be called from bash (as
% opposed to the wrapped function).
%
%
% ARGUMENTS
% =========
% NOTE: See implementation and bicas.tools.batch.main() for details.
% --
% varargin : Paths to input datasets, or directories with datasets (recursive).
%       input paths
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-03-19.
%
function main_bash(bicasConfigFile, isCdagOption, modeStr, outputDir, referenceDir, varargin)

switch(isCdagOption)
  case '-n'    % N = normal/no CDAG
    outputIsCdag = false;
  case '-c'    % C = CDAG
    outputIsCdag = true;
  otherwise
    error('Illegal isCdagOption="%s"', isCdagOption)
end

bicas.tools.batch.main(...
  bicasConfigFile, outputIsCdag, modeStr, outputDir, ...
  referenceDir, varargin)
end
