%
% Implementation for being able to build SWMs for tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef TestSwmProcessing < bicas.proc.SwmProcessing



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
    end    % properties(SetAccess=immutable)



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        % ARGUMENTS
        % =========
        % InputsMap
        %       containers.Map with
        %       <keys>   : String defining a name of an input ("prodFuncInputKey" in
        %                  bicas.swm.SoftwareModeList).
        %       <values> : A struct with data corresponding to a CDF file
        %                  (zVariables+global attributes).
        % OutputsMap
        %       containers.Map with
        %       <keys>   : String defining a name of an output ("prodFuncOutputKey" in
        %                  bicas.swm.SoftwareModeList).
        %       <values> : A struct with data corresponding to a CDF file (zVariables).
        %
        % OVERRIDE
        function OutputDatasetsMap = production_function(obj, ...
            InputDatasetsMap, rctDir, NsoTable, Bso, L)

            OutputDatasetsMap = containers.Map();
        end



    end    % methods(Access=public)



end
