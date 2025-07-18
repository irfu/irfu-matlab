%
% Class that represents how to convert one particular QRCID into modifications
% of quality ZVs.
%
% NOTE: The class does not include the QRCID itself.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSetting
  % PROPOSAL: Better name
  %   setting
  %       CON: Could be confused with more proper settings such as the NSO
  %            table itself. The information stored in the object is more like
  %            a constant.
  %   interpretation
  %   configuration
  %   behaviour
  %   entry
  %   policy
  %   entry
  %   action
  %   QRC, QRCID
  %   quality ZVs
  %   quality ZVs modification
  %   --
  %   QRCC = QRC Configuration
  %   QRCS = QRC Setting -- IMPLEMENTED
  %   --
  %   PROPOSAL: Abbreviation
  %
  % PROPOSAL: Include UFV as an action.
  %   PRO: Can redefine excluding data using BDM using this class.
  %     CON: Name Quality-Related Condition becomes a misnomer.
  %     CON: Should probably phase out that functionality anyway.
  %   PRO: Useful for blanking antenna fails.
  %     https://github.com/irfu/irfu-matlab/issues/142
  %   CON/PROBLEM: Needs way of specifying channel. All channels?!!
  %     CON-PROPOSAL: Include "channel addresses"?
  %
  % PROPOSAL: Modify definition of QRC beyond quality-related conditions to
  %           include any special treatment of data.
  %   Ex: BIAS current sweeps.
  %   Ex: BIAS off (BW).
  %   See readme.txt on special cases.
  %   TODO-DEC: Term to replace "QRC".
  %     ~Special case
  %     ~Special treatment/processing/handling
  %     ~Data treatment/processing/handling
  %     ~Exception
  %       CON: May imply programming term exception=error handling.
  %     ~Condition
  %     Special Case Condition
  %     Special Condition
  %
  % PROPOSAL: Use only keyword arguments with default values.
  %   PRO: Common to not set all values.
  %     Ex: Hardcoding tests.
  %   CON: Longer calls.



  properties(SetAccess=immutable)
    % *Cap* (max value) for the CDF ZV "QUALITY_FLAG" when the QRC applies.
    %
    % NOTE: This is interpretation is in compliance with how the ZV
    % QUALITY_FLAG is supposed to be set/updated.
    QUALITY_FLAG

    % Bits that should be set in either ZV "L2_QUALITY_BITMASK" or
    % "L3_QUALITY_BITMASK". The context in which this class is used
    % determines which.
    % NOTE: The value is supposed to be OR:ed with a preceding value, i.e. only
    % set bits override the previous value.
    Lx_QUALITY_BITMASK
  end



  methods(Access=public)



    function obj = QrcSetting(QUALITY_FLAG, Lx_QUALITY_BITMASK)

      assert(bicas.utils.validate_ZV_QUALITY_FLAG(QUALITY_FLAG))
      assert(isa(Lx_QUALITY_BITMASK, 'uint16'))

      obj.QUALITY_FLAG       = QUALITY_FLAG;
      obj.Lx_QUALITY_BITMASK = Lx_QUALITY_BITMASK;
    end



  end



end
