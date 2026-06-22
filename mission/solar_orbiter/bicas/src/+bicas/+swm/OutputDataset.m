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
    dsid
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
        cliOptionHeaderBody, dsid, pfoid, ...
        swdName, swdDescription, skeletonVersion)

      [~, datasetLevel, ~] = solo.adm.disassemble_DATASET_ID(dsid);

      obj.cliOptionHeaderBody = cliOptionHeaderBody;
      obj.dsid                = dsid;

      obj.pfoid               = pfoid;
      obj.swdName             = swdName;
      obj.swdDescription      = swdDescription;
      obj.skeletonVersion     = skeletonVersion;

      % NOTE: Not argument. Extracted from DSID.
      obj.datasetLevel        = datasetLevel;

      bicas.swm.utils.assert_SIP_CLI_option(obj.cliOptionHeaderBody)
      bicas.swm.utils.assert_text(          obj.swdName)
      bicas.swm.utils.assert_text(          obj.swdDescription)
      bicas.swm.utils.assert_DSID(          obj.dsid)
      solo.adm.assert_dataset_level(        obj.datasetLevel)
      bicas.assert_skeleton_version(        obj.skeletonVersion)
    end



  end    % methods(Access=public)



end
