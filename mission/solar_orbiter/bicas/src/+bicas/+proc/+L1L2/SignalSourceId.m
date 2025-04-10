%
% Immutable class which instances represent the source of a signal, i.e.
% either:
% (1) an ASR (ASID), or
% (2) various special cases.
%
% NOTE: Can not represent the source of a reconstructed signal, e.g. a diff
% calculated by subtracting to (calibrated) singles.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SignalSourceId
  % PROPOSAL: Use for solo.BSACT_utils.
  %
  % PROPOSAL: Should include three separate cases for 2.5V_REF and GND
  %           respectively (i.e. 2 cases --> 2x3 cases).



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=public)
    Asid
  end
  properties(SetAccess=immutable, GetAccess=private)
    % NOTE: Private value. Value can (and should) be indirectly accessed by
    %       comparing the object (isequaln) with one of the object constants.
    specialCase
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
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



end
