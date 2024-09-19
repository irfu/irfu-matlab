%
% Immutable class which instances represent the destination of a signal, i.e.
% where the data should ultimately be stored in the output datasets, i.e.
% either:
% (1) an ASR (to determine which zVariable to store it in in the output
%     dataset), or
% (2) "nowhere" (when BDM is unknown).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SignalDestinationId
  % PROPOSAL: Better name.
  %   (output) dataset, zVariables



  properties(GetAccess=public, Constant)
    C = bicas.proc.L1L2.SignalDestinationId.init_const();
  end



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=public)
    % ASID object or empty.
    Asid

    % Whether destination is "NOWHERE", i.e. the signal does not have a
    % destination and should be ignored.
    isNowhere
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Constructor
    function obj = SignalDestinationId(value)
      if isa(value, 'bicas.proc.L1L2.AntennaSignalId')
        obj.Asid      = value;
        obj.isNowhere = false;
      elseif isequal(value, 'NOWHERE')
        obj.Asid      = [];
        obj.isNowhere = true;
      else
        error('BICAS:Assertion:IllegalArgument', 'Illegal argument.')
      end
    end



  end    % methods(Access=public)



  methods(Access=private, Static)



    function C = init_const()
      C = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants( ...
        @(Asid) (bicas.proc.L1L2.SignalDestinationId(Asid)));

      C.NOWHERE = bicas.proc.L1L2.SignalDestinationId('NOWHERE');
    end



  end



end
