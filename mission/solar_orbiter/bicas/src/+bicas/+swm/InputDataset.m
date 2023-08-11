%
% Class that stores metadata for one SWD input dataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef InputDataset

    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        cliOptionHeaderBody
        prodFuncInputKey
        datasetId
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)

        function obj = InputDataset(...
                cliOptionHeaderBody, datasetId, prodFuncInputKey)

            % NOTE: No dataset/skeleton version.
            obj.cliOptionHeaderBody = cliOptionHeaderBody;
            obj.prodFuncInputKey    = prodFuncInputKey;
            obj.datasetId           = datasetId;

            bicas.swm.utils.assert_SIP_CLI_option(obj.cliOptionHeaderBody)
            % NOTE: Using the INTERNAL assertion function, not the global one.
            bicas.swm.utils.assert_DATASET_ID(    obj.datasetId)
        end

    end

end
