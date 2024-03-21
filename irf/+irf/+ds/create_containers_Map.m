%
% Create a containers.Map object in a way that works around MATLAB shortcomings.
% Specifically:
% -- Initialize containers.Map and simultaneosly setting keys+values and key+value types.
% -- Initialize containers.Map with empty key+value arrays (special case behaves differently in constructor).
%
%
% NOTE: Forces user to specify key+value type.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-07-23.
%
function Map = create_containers_Map(keyType, valueType, keyCa, valueCa)
Map = containers.Map('KeyType', keyType, 'ValueType', valueType);

assert(iscell(keyCa))
assert(iscell(valueCa))
assert(numel(keyCa) == numel(valueCa), ...
  'Arguments keyCa and valueCa do not have same length.')

for i = 1:numel(keyCa)
  Map(keyCa{i}) = valueCa{i};
end
end
