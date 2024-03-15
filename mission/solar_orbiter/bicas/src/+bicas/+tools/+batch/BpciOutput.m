%
% Helper class to bicas.tools.batch.BicasProcessingCallInfo.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef BpciOutput
  % PROPOSAL: Replace property "filename" with "path". -- IMPLEMENTED
  %   PRO: Analogous with bicas.tools.batch.BpciInput.
  %   PRO: Can maybe use object in BPTD (instead of this class converted to
  %        struct which is then modified).
  %   PRO: A path is actually sent to BICAS. More accurate representation of
  %        a call to BICAS.
  %   NOTE: Autocreation of BPCIs needs way of setting paths (not just
  %         filenames) of output datasets.
  %       Ex: bicas.tools.batch.autocreate_one_SWM_BPCI().
  %       PROPOSAL: Submit output directory to those functions.
  %       PROPOSAL: Redefine ~createOutputFilenameFh to set entire path.
  %           PRO: Most generic.
  %   NOTE/PROBLEM: Class becomes identical to
  %                 bicas.tools.batch.BpciInput!



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    cohb
    dsi
    path
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = BpciOutput(cohb, dsi, path)
      obj.cohb = cohb;
      obj.dsi  = dsi;
      obj.path = path;
    end



  end    % methods(Access=public)



end
