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
    dsid
    pfiid
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = InputDataset(cliOptionHeaderBody, dsid, pfiid)

      % NOTE: No dataset/skeleton version.
      obj.cliOptionHeaderBody = cliOptionHeaderBody;
      obj.dsid                 = dsid;
      obj.pfiid               = pfiid;

      bicas.swm.utils.assert_SIP_CLI_option(obj.cliOptionHeaderBody)
      bicas.swm.utils.assert_DSID(           obj.dsid)
    end



  end



end
