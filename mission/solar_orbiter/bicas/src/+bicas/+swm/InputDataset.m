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
        dsi
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)

        function obj = InputDataset(...
                cliOptionHeaderBody, dsi, prodFuncInputKey)

            % NOTE: No dataset/skeleton version.
            obj.cliOptionHeaderBody = cliOptionHeaderBody;
            obj.prodFuncInputKey    = prodFuncInputKey;
            obj.dsi           = dsi;

            bicas.swm.utils.assert_SIP_CLI_option(obj.cliOptionHeaderBody)
            % NOTE: Using the INTERNAL assertion function, not the global one.
            bicas.swm.utils.assert_DSI(    obj.dsi)
        end

    end

end
