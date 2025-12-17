%
% Helper class to bicas.tools.batch.BicasProcessingCallInfo.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef BpciInput



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    cohb
    dsid
    path
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = BpciInput(cohb, dsid, path)
      obj.cohb = cohb;
      obj.dsid  = dsid;
      obj.path = path;
    end



  end    % methods(Access=public)



end
