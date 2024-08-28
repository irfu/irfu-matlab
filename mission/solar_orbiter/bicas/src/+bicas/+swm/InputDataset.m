%
% Class that stores metadata for one SWD input dataset.
%
% NOTE: Not to be confused with bicas.InputDataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef InputDataset
  % PROPOSAL: Rename to not be same as bicas.InputDataset and bicas.OutputDataset.
  %   InputSwmDataset
  %       CON: SWM already in path (MATLAB package name).
  %   InputMetadataDataset, InputDatasetMetadata
  %       CON: 2x "data"
  % PROPOSAL: Use abbreviations which are different from bicas.InputDataset
  %           abbrevations.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    cliOptionHeaderBody
    pfiid
    dsi
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = InputDataset(...
        cliOptionHeaderBody, dsi, pfiid)

      % NOTE: No dataset/skeleton version.
      obj.cliOptionHeaderBody = cliOptionHeaderBody;
      obj.pfiid               = pfiid;
      obj.dsi                 = dsi;

      bicas.swm.utils.assert_SIP_CLI_option(obj.cliOptionHeaderBody)
      bicas.swm.utils.assert_DSI(           obj.dsi)
    end



  end



end
