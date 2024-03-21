%
% Given two N-dimensional arrays (y1, y2), with associated coordinate values in
% each dimension (x1Ca, x2Ca), create a third N-dimensional array (Y) with
% associated coordinate values in each dimension (xCa). The return values are
% the union of the arguments, with a fill value used for unset values.
%
% NOTE: "coordinated" in the function name refers to N-dimensional arrays where
% one (Cartesian) coordinate axis is associated with every array dimension, and
% one coordinate axis value is associated with every index value. Therefore
% every array component also has a coordinate.
%
%
% ARGUMENTS
% =========
% fillVal    : Scalar numeric. Value used for components in y not determined by
%              y1 and y2.
% x1Ca, x2Ca : Size Nx1 cell arrays of 1D numeric arrays. {i} must not be empty.
% y1, y2     : Numeric arrays of size length(x1Ca{1}) x ... x length(x1Ca{N}).
%              and                    length(x2Ca{1}) x ... x length(x2Ca{N}).
%
%
% RETURN VALUES
% =============
% xCa : Column cell array. xCa{i} = union(x1Ca{i}, x2Ca{i}).
% y   : Numeric array of size length(xCa{1}) x ... x length(xCa{N}).
%       Values come from y1, and y2, and fillVal. In case of collisions, then y2
%       takes precedence over y1.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-08-14.
%
function [xCa, y] = merge_coordinated_arrays(fillVal, x1Ca, y1, x2Ca, y2)
% PROPOSAL: Better name?
% NOTE: Does not guarantee/warn of overlap between y1 and y2. y1 values could be overwritten by y2.

% IMPLEMENTATION NOTE: Every N-dimensional array can in principle be seen as
% an infinite-dimensional array, with size=1 for dimensions >=N+1. These
% higher dimensions do not have coordinates associated with them. A
% zero-dimensional array is thus equivalent to a scalar, without any
% coordinates, and it is distinct from e.g. a 2D array of size 1x1, which
% would have coordinate values for the first and second dimension.
% Similarily, an empty array, e.g. 0x0, therefore has at least dimension.

% ASSERTIONS
assert(isscalar(fillVal))
nDims = irf.assert.sizes(x1Ca(:), [-1], x2Ca(:), [-1]);



% NOTE: Used as argument in ones(ySize). Must therefore be
%   (1) row vector, and
%   (2) at least size 1x2 (with default value 1).
ySize = ones(1, max(nDims, 2));

xCa   = cell(nDims, 1);
S1.type = '()';
S2.type = '()';
S1.subs = {};
S2.subs = {};
for i = 1:nDims
  % ASSERTIONS
  % There is no known fundamental reason to not permit zero coordinates on
  % axis. The code is just not vetted for that case yet.
  assert(length(x1Ca{i}) >= 1)
  assert(length(x2Ca{i}) >= 1)
  assert(length(x1Ca{i}) == size(y1, i), 'Size of y1 in dimension %i does not match size of x1Ca{%i}.', i, i)
  assert(length(x2Ca{i}) == size(y2, i), 'Size of y2 in dimension %i does not match size of x2Ca{%i}.', i, i)
  % Coordinates must be unique for assignment to be unique.
  irf.assert.number_set(x1Ca{i})
  irf.assert.number_set(x2Ca{i})



  xCa{i} = union(x1Ca{i}, x2Ca{i});    % NOTE: Result is sorted.
  xCa{i} = xCa{i}(:);    % Force column array (for testing).
  ySize(i) = length(xCa{i});

  [~, j1] = ismember(x1Ca{i}, xCa{i});
  [~, j2] = ismember(x2Ca{i}, xCa{i});

  S1.subs{i} = j1;
  S2.subs{i} = j2;
end

y = ones(ySize) * fillVal;

% NOTE: Empirically, subsasgn can reliably crash MATLAB when y=NaN, y1=[],
% and S1.subs={}, i.e. zero-dimensions
% !!! Must therefore treat this as a special case.
if nDims == 0
  y = y2;
else
  y = subsasgn(y, S1, y1);
  y = subsasgn(y, S2, y2);
end
end
