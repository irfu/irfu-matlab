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
function Map = create_containers_Map(keyType, valueType, keyList, valueList)
    Map = containers.Map('KeyType', keyType, 'ValueType', valueType);
    
    assert(iscell(keyList))
    assert(iscell(valueList))
    assert(numel(keyList) == numel(valueList))
    
    for i = 1:numel(keyList)
        Map(keyList{i}) = valueList{i};
    end
end
