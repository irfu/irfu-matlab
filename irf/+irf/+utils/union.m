%
% Set union, as defined by MATLAB function "union" (default arguments), but
% applied to an arbitrary number of sets.
%
%
% ARGUMENTS
% =========
% varargin: n>=1 arguments representing sets.
%
%
% RETURN VALUES
% =============
% u : The union of all the input datasets as defined by MATLAB function "union"
%     1D column array. Only unique values. Sorted.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-09-24.
%
function u = union(varargin)
nSets = numel(varargin);

% IMPLEMENTATION NOTE: One could include case nSets==1 in nSets>=1 by
% starting with u=empty set, but one must then create an empty set whose
% type (cell, numeric, logical etc) depending on the other sets which makes
% the implementation uglier than a special case.

if nSets == 1
  % IMPLEMENTATION NOTE: Special treatment to be consistent with the
  % output of "union".
  % NOTE: Always 1D (not unique in every row och column etc).
  u = unique(varargin{1});   % Will automatically be sorted.
elseif nSets >= 2
  u = varargin{1};
  for i = 2:nSets
    % NOTE: Always 1D (not union of every row och column etc).
    u = union(u, varargin{i});
  end
else
  % IMPLEMENTATION NOTE: In principle, zero input sets should yield an empty
  % output set, but it ambiguous what kind of set: cell array, numeric array,
  % logical, char?
  error('No input sets. Function is not well-defined for this case.')
end

% For consistency.
u = u(:);
end
