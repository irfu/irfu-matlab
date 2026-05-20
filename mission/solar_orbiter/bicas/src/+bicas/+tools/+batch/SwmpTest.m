%
% Implementation for being able to build SWMs for tests.
%
% NOTE: The class is not configurable and does not do anything.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SwmpTest < bicas.proc.SwmProcessing



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % OVERRIDE
    function OutputDatasetsMap = production_function(obj, ...
        InputDatasetsMap, rctDir, NsoTable, Bso, L)

      OutputDatasetsMap = containers.Map();
    end



  end    % methods(Access=public)



end
