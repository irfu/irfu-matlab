%
% Wrapper around irf.utils.split_by_change() to make the return value on the
% (almost) same format as bicas.utils.group_unique_rows() to make it easy to
% switch between the two.
%
%
% ARGUMENTS
% =========
% varargin
%       See irf.utils.split_by_change().
%
%
% RETURN VALUES
% =============
% iGroupArCa
%       Column cell array of row indices.
%       {iGroup} = Sorted column array of indices to rows in arguments.
%       NOTE: The sequence of values {iGroup}(1) is sorted.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function iGroupArCa = group_by_change(varargin)

[i1Ar, i2Ar, n] = irf.utils.split_by_change(varargin{:});

iGroupArCa = cell(n, 1);
for i = 1:n
  iGroupArCa{i, 1} = [i1Ar(i):i2Ar(i)]';
end

end
