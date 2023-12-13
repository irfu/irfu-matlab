%
% Class that stores metadata and production function for one SWM.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SoftwareMode
    % PROPOSAL: Fieldname change
    %   inputsList  --> inputsArray
    %   outputsList --> outputsArray
    %   NOTE: Likely influences BICAS testing code and pipeline.
    %         Should only be implemented at the right time.

    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        prodFunc
        cliOption
        swdPurpose
        inputsList
        outputsList
    end    % properties(SetAccess=immutable)



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = SoftwareMode(...
                prodFunc, cliOption, swdPurpose, ...
                inputsList, outputsList)

            obj.prodFunc    = prodFunc;
            % NOTE: s/w mode CLI _ARGUMENT_ is not intended to be prefixed by
            % e.g. "--". Variable therefore NOT named *Body.
            obj.cliOption   = cliOption;
            obj.swdPurpose  = swdPurpose;
            obj.inputsList  = inputsList;
            obj.outputsList = outputsList;

            %============
            % ASSERTIONS
            %============
            bicas.swm.utils.assert_SWM_CLI_option(obj.cliOption)
            bicas.swm.utils.assert_text(          obj.swdPurpose)

            % Important. Check uniqueness of SIP options.
            irf.assert.castring_set( {...
                obj.inputsList( :).cliOptionHeaderBody, ...
                obj.outputsList(:).cliOptionHeaderBody })

            assert(isa(obj.inputsList,  'bicas.swm.InputDataset'))
            assert(isa(obj.outputsList, 'bicas.swm.OutputDataset'))

            irf.assert.castring_set( { obj.inputsList(:).prodFuncInputKey   })
            irf.assert.castring_set( { obj.outputsList(:).prodFuncOutputKey })
        end



    end    % methods(Access=public)

end
