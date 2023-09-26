%
% Dumb class for making it possible to have a handle to a non-handle object,
% e.g. regular arrays. Can be used for ensuring meaningful pre-allocation inside
% e.g. containers.Map values to increase performance.
%
% Used by bicas.utils.SameRowsMap for increasing performance.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef HandleWrapper < handle
    
    properties(Access=public)
        v
    end
    
    methods
        function obj = HandleWrapper(v)
            obj.v = v;
        end
    end
end