%
% Call strrep multiple times for the same string.
%
%
% ARGUMENTS
% =========
% varargin
%       Pairs of strings: old substring to search for and its replacement
%       substring.
%       NOTE: The replacement substring may in turn be replaced by another
%       substring. Therefore, the order of specified replacements matters.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-12-03
%
function s = strrep_many(s, varargin)
% PROPOSAL: Change argument syntax (varargin) to something that is more suitable for data structures AND hardcoding.
%   PROPOSAL: oldSsList, newSsList
%       CON: Bad for hardcoding. Each string pair can not be on the same row.
%   PROPOSAL: replacementTable: {i}{1} = old substring, {i}{2} = new substring
%   PROPOSAL: replacementTable: {i, 1} = old substring, {i, 2} = new substring
%       PRO: oldSsList, newSsList can easily be merged to such a 2D cell array.

while true
  if numel(varargin) >= 2
    oldSs = varargin{1};
    newSs = varargin{2};

    % ASSERTIONS
    irf.assert.castring(oldSs)
    irf.assert.castring(newSs)

    varargin = varargin(3:end);
  elseif numel(varargin) == 0
    return
  else
    error('BICAS:Assertion', ...
      ['Not same number of arguments for old substrings and', ...
      ' new substrings.'])
  end

  s = strrep(s, oldSs, newSs);
end
end
