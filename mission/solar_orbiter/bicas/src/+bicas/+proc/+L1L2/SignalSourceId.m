%
% Immutable class which instances represent the source of a signal, i.e.
% either:
% (1) an ASR (ASID), or
% (2) various special cases.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SignalSourceId
  % PROPOSAL: Use for solo.BSACT_utils.



  properties(GetAccess=public, Constant)
    C = bicas.proc.L1L2.SignalSourceId.init_const()
  end



  properties(SetAccess=immutable, GetAccess=public)
    Asid
  end



  properties(SetAccess=immutable, GetAccess=private)
    % NOTE: Private value. Value can (and should) be indirectly accessed by
    %       comparing the object (isequaln) with one of the object constants.
    specialCase
  end



  methods(Access=public)



    % Constructor
    function obj = SignalSourceId(value)
      if isa(value, 'bicas.proc.L1L2.AntennaSignalId')
        obj.Asid        = value;
        obj.specialCase = [];
      elseif ischar(value) && ismember(value, {'2.5V_REF', 'GND', 'UNKNOWN'})
        obj.Asid        = [];
        obj.specialCase = value;
      else
        error('BICAS:Assertion:IllegalArgument', 'Illegal argument.')
      end
    end



    function isAsr = is_ASR(obj)
      isAsr = isa(obj.Asid, 'bicas.proc.L1L2.AntennaSignalId');
    end



  end    % methods(Access=public)



  methods(Access=private, Static)



    function C = init_const()
      C = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants( ...
        @(Asid) (bicas.proc.L1L2.SignalSourceId(Asid)));

      C.REF25V   = bicas.proc.L1L2.SignalSourceId('2.5V_REF');
      C.GND      = bicas.proc.L1L2.SignalSourceId('GND');
      C.UNKNOWN  = bicas.proc.L1L2.SignalSourceId('UNKNOWN');
    end



  end



end
