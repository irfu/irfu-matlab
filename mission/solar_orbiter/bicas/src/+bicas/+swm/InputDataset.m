%
% Class that stores metadata for one SWD input dataset.
% 
% NOTE: Not to be confused with bicas.InputDataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef InputDataset
% PROPOSAL: Rename to not be same as bicas.InputDataset.
%   InputSwmDataset
%       CON: SWM already in path (MATLAB package name).
%   InputMetadataDataset, InputDatasetMetadata
%       CON: 2x "data"



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
            obj.dsi                 = dsi;

            bicas.swm.utils.assert_SIP_CLI_option(obj.cliOptionHeaderBody)
            bicas.swm.utils.assert_DSI(           obj.dsi)
        end



    end



end
