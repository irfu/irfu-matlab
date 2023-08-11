%
% Miscellaneous utility functions in the form of static methods.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils

    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)

        % Assert that string contains human-readable text.
        function assert_text(str)
            irf.assert.castring_regexp(str, '.* .*')
            irf.assert.castring_regexp(str, '[^<>]*')
        end



        function assert_SWM_CLI_option(swmCliOption)
            irf.assert.castring_regexp(...
                swmCliOption, ...
                bicas.constants.SWM_CLI_OPTION_REGEX)
        end



        % NOTE: Really refers to "option body".
        function assert_SIP_CLI_option(sipCliOptionBody)
            irf.assert.castring_regexp(...
                sipCliOptionBody, ...
                bicas.constants.SIP_CLI_OPTION_BODY_REGEX)
        end



        % NOTE: Wrapper around global counterpart.
        function assert_DATASET_ID(datasetId)
            bicas.assert_BICAS_DATASET_ID(datasetId)

            % ASSERTION: Only using SOLO_* DATASET_IDs.
            [sourceName, ~, ~] = solo.adm.disassemble_DATASET_ID(datasetId);
            assert(strcmp(sourceName, 'SOLO'))
        end

    end    % methods(Static)

end
