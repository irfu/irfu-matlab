%
% Class that stores metadata for one SWD output dataset.
%
% NOTE: Not to be confused with bicas.OutputDataset
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef OutputDataset
  % PROPOSAL: Rename
  %   See bicas.swm.InputDataset comments.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    cliOptionHeaderBody
    dsi
    datasetLevel
    pfoid

    % "name" in SWD: Human-readable name, but shorter than swdDescription.
    swdName
    % "description in SWD: Human-readable description.
    swdDescription

    skeletonVersion
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = OutputDataset(...
        cliOptionHeaderBody, dsi, pfoid, ...
        swdName, swdDescription, skeletonVersion)

      [~, datasetLevel, ~] = solo.adm.disassemble_DATASET_ID(dsi);

      obj.cliOptionHeaderBody = cliOptionHeaderBody;
      obj.dsi                 = dsi;
      obj.datasetLevel        = datasetLevel;

      obj.pfoid               = pfoid;
      obj.swdName             = swdName;
      obj.swdDescription      = swdDescription;
      obj.skeletonVersion     = skeletonVersion;

      bicas.swm.utils.assert_SIP_CLI_option(obj.cliOptionHeaderBody)
      bicas.swm.utils.assert_text(          obj.swdName)
      bicas.swm.utils.assert_text(          obj.swdDescription)
      bicas.swm.utils.assert_DSI(           obj.dsi)
      solo.adm.assert_dataset_level(        obj.datasetLevel)
      bicas.assert_skeleton_version(        obj.skeletonVersion)
    end



  end    % methods(Access=public)



end
