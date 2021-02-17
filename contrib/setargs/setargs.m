function argstruct = setargs(defaultargs, varargs)
% SETARGS Name/value parsing and assignment of varargin with default values
% 
% This is a utility for setting the value of optional arguments to a
% function. The first argument is required and should be a cell array of
% "name, default value" pairs for all optional arguments. The second
% argument is optional and should be a cell array of "name, custom value"
% pairs for at least one of the optional arguments.
% 
%   USAGE: argstruct = setargs(defaultargs, varargs)
% __________________________________________________________________________
% OUTPUT
% 
% 	ARGSTRUCT
%    structure containing the final argument values
% __________________________________________________________________________
% INPUTS
% 
% 	DEFAULTARGS  
%     cell array of "'Name', value" pairs for all variables with default
%     values
% 
% 	VARARGS [optional]     
%     cell array of user-specified "'Name', value" pairs for one or more of
%     the variables with default values. this will typically be the
%     "varargin" cell array. for each pair, SETARGS determines if the
%     specified variable name can be uniquely matched to one of the default
%     variable names specified in DEFAULTARGS. matching uses STRNCMPI and
%     thus is case-insensitive and open to partial name matches (e.g.,
%     default variable name 'FontWeight' would be matched by 'fontweight',
%     'Fontw', etc.). if a match is found, the user-specified value is then
%     used in place of the default value. if no match is found or if
%     multiple matches are found, SETARGS returns an error and displays in
%     the command window information about the argument that caused the
%     problem.
% __________________________________________________________________________
% USAGE EXAMPLE (TO BE USED AT TOP OF FUNCTION WITH VARARGIN)
% 
%     defaultargs = {'arg1', 0, 'arg2', 'words', 'arg3', rand}; 
%     argstruct   = setargs(defaultargs, varargin)
%


% ---------------------- Copyright (C) 2015 Bob Spunt -----------------------
%	Created:  2015-03-11
%	Email:    spunt@caltech.edu
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or (at
%   your option) any later version.
%       This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%       You should have received a copy of the GNU General Public License
%   along with this program.  If not, see: http://www.gnu.org/licenses/.
% 
if nargin < 1, mfile_showhelp; return; end
if nargin < 2, varargs = []; end
defaultargs = reshape(defaultargs, 2, length(defaultargs)/2)'; 
if ~isempty(varargs)
    if mod(length(varargs), 2)
        error('Optional inputs must be entered as "''Name'', Value" pairs, e.g., myfunction(''arg1'', val1, ''arg2'', val2)'); 
    end
    arg = reshape(varargs, 2, length(varargs)/2)';
    for i = 1:size(arg,1)
       idx = strncmpi(defaultargs(:,1), arg{i,1}, length(arg{i,1}));
       if sum(idx) > 1
           error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaultargs{idx, 1});
       elseif ~any(idx)
           error('Input "%s" does not match a valid input.', arg{i,1});
       else
           defaultargs{idx,2} = arg{i,2};
       end
    end
end

for i = 1:size(defaultargs,1), assignin('caller', defaultargs{i,1}, defaultargs{i,2}); end
if nargout>0, argstruct = cell2struct(defaultargs(:,2), defaultargs(:,1)); end
end
% =========================================================================
% * SUBFUNCTIONS
% =========================================================================
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));  
end