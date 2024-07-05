%
% Represents a concrete output dataset produced by processing, primarily the
% content.
%
% NOTE: Not to be confused with bicas.swm.OutputDataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef OutputDataset
  % PROPOSAL: Automatic test code.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    Zv
    Ga
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = OutputDataset(Zv, Ga)
      obj.Zv = Zv;
      obj.Ga = Ga;
    end



  end    % methods(Access=public)



end
