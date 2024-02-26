%
% Helper class to bicas.tools.batch.BicasProcessingCallInfo.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef BpciInput



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        cohb
        dsi
        path
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = BpciInput(cohb, dsi, path)
            obj.cohb = cohb;
            obj.dsi  = dsi;
            obj.path = path;
        end



    end    % methods(Access=public)



end
