%
% Help function or initializing an empty struct. Useful for initializing an
% empty struct array with fieldnames, and of a given size (e.g. 0x1 instead of
% 0x0).
%
%
% ARGUMENTS
% =========
% size
%       1D numeric array. If size implies a non-zero number of elemnts, then the
%       fields will all be set to [].
% varargin
%       Arbitrary number of arguments representing fieldnames.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-17.
%
function S = empty_struct(size, varargin)
%
% PROPOSAL: Be able to set the field value (same for all fields) for non-empty stuct arrays.
%   Ex: {}, cell(0,1)

% NOTE: cell(n) ==> nxn cell array.
assert(isnumeric(size), 'Argument "size" is not numeric.')
irf.assert.vector(size)
size(end+1) = 1;
size(end+1) = 1;
%assert(numel(size) >= 2, '')

if isempty(varargin)
  % CASE: Create struct with no fieldnames

  % IMPLEMENTATION NOTE: Not obvious how to create an EMPTY struct ARRAY
  % without fields. Uses trick to achieve it.
  S = rmfield(struct('a', cell(size)), 'a');
else
  % CASE: Create struct with fieldnames

  structArgs = {};
  for i = 1:numel(varargin)
    structArgs{2*i-1} = varargin{i};
    structArgs{2*i  } = cell(size);
  end

  S = struct(structArgs{:});
end
end
