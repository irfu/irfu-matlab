%
% Class that stores metadata and production function for one SWM.
%
% NOTE: Not to be confused with class bicas.proc.SwmProcessing which contains
% the processing associated with a SWM but not the metadata for the SWD. This
% class does however contain a copy to the relevant instance of
% bicas.proc.SwmProcessing.
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
    Swmp
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
        Swmp, cliOption, swdPurpose, ...
        inputsList, outputsList)

      obj.Swmp        = Swmp;
      % NOTE: s/w mode CLI _ARGUMENT_ is not intended to be prefixed by
      % e.g. "--". Variable is therefore NOT named *Body.
      obj.cliOption   = cliOption;
      obj.swdPurpose  = swdPurpose;
      obj.inputsList  = inputsList;
      obj.outputsList = outputsList;

      %============
      % ASSERTIONS
      %============
      assert(isa(obj.Swmp, 'bicas.proc.SwmProcessing'))
      bicas.swm.utils.assert_SWM_CLI_option(obj.cliOption)
      bicas.swm.utils.assert_text(          obj.swdPurpose)

      % Important. Check uniqueness of SIP options.
      irf.assert.castring_set( {...
        obj.inputsList( :).cliOptionHeaderBody, ...
        obj.outputsList(:).cliOptionHeaderBody })

      assert(isa(obj.inputsList,  'bicas.swm.InputDataset'))
      assert(isa(obj.outputsList, 'bicas.swm.OutputDataset'))

      irf.assert.castring_set( { obj.inputsList(:).pfiid   })
      irf.assert.castring_set( { obj.outputsList(:).pfoid })
    end



  end    % methods(Access=public)

end
