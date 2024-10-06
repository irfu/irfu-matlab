%
% Immutable class which instances represent an SSID with metadata.
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
  % PROPOSAL: Should include three separate cases for REF25V and GND
  %           respectively (i.e. 2 cases --> 2x3 cases).



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=public)
    % ASID or empty.
    asid
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
    function obj = SignalSourceId(asidOrSpecialCase)
      if isa(asidOrSpecialCase, 'uint8')
        obj.asid        = asidOrSpecialCase;
        obj.specialCase = [];
      elseif isstring(asidOrSpecialCase) ...
          && ismember(asidOrSpecialCase, ["REF25V", "GND", "UNKNOWN"])
        obj.asid        = [];
        obj.specialCase = asidOrSpecialCase;
      else
        error('BICAS:Assertion:IllegalArgument', 'Illegal argument.')
      end
    end



    function isAsr = is_ASR(obj)
      isAsr = ~isempty(obj.asid);
    end



  end    % methods(Access=public)



end
