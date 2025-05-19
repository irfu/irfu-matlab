%
% Handle class
%
% Class which effectively wraps a dictionary
% (ASR) SDID-->bicas.proc.L1L2.SdChannelData.
%
% This useful since
% (1) (MATLAB) dictionary values can not be non-scalar
%     (bicas.proc.L1L2.SdChannelData is a column vector). Therefore,
%     the implementation must work around this (it uses 1x1 cell arrays).
% (2) it can sum up the number of SDCD fill positions.
%
% NOTE: The constructor does not initialize the object completely because
% the constructor call would be too large and awkward then.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SdChannelDataDict < handle
  % PROPOSAL: Better name
  %   ~ASR (as opposed to non-ASR)
  %   ~9x ASR (as opposed to 5x BLTS)
  %   ~9x ASR (as opposed to just one channel)
  %   ~SdChannelData
  %     Which may be used for BLTS(?)
  %   ~Dictionary
  %   ~Collection
  %   --
  %   ASDCD = AsrSdChannelData
  %     CON: Does not imply collection of ALL ASR channels.
  %
  % PROPOSAL: Better tests.
  %
  % PROPOSAL: Constructor pre-allocates SDCDs.
  % PROPOSAL: Implement custom print version of class.
  %   PRO: Useful for debugging.
  %
  % PROPOSAL: Constructor must set values for all keys.
  %   PRO: Class is always initialized.
  %   PROPOSAL: Cell array, 9x2, (iSsid, 1)=SSID, (iSsid, 2)=SDCD
  %     CON: Memory problems?
  %   PRO: More natural to assert similarities between constituent SDCDs.
  %
  % PROPOSAL: Assert that all SDCDs have the same size and data types.



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant, Access=private)
    % ASR SDIDs.
    PERMITTED_KEYS_AR = bicas.proc.L1L2.const.C.SDID_ASR_AR
  end
  properties(Dependent)
    % Total number of bWholeRowIsNan values in all the SDCD objects stored
    % within this object combined.
    nWholeRowIsNan
  end



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Access=private)
    Dict
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = SdChannelDataDict()
      obj.Dict = dictionary;
    end



    % Set key value.
    function obj = set(obj, asrSdid, Sdcd)
      assert(isscalar(asrSdid))
      assert(ismember(asrSdid, obj.PERMITTED_KEYS_AR))
      assert(isa(Sdcd, 'bicas.proc.L1L2.SdChannelData'))

      obj.Dict(asrSdid) = {Sdcd};
    end



    function Sdcd = get(obj, asrSdid)
      assert(isscalar(asrSdid))
      ca   = obj.Dict(asrSdid);
      Sdcd = ca{1};
    end



    % function groups = getPropertyGroups(obj)
    %   % PROPOSAL: Separate properties for MATLAB class and size.
    %   %   PRO: Avoids repetition.
    %   %   CON: Less good for debugging class itself.
    %
    %   % IMPLEMENTATION NOTE: It appear that one can only represent
    %   % "properties" using single-row strings.
    %
    %   properties = struct(...
    %     'dataAr', bicas.utils.FPArray.value_to_single_row_string(obj.dataAr, obj.fpAr), ...
    %     'fpAr',   bicas.utils.FPArray.value_to_single_row_string(obj.fpAr), ...
    %     'size',   size(obj), ...
    %     'mc',     obj.mc, ...
    %     'onlyFp', all( obj.fpAr, 'all'), ...
    %     'noFp',   all(~obj.fpAr, 'all') ...
    %     );
    %   groups = matlab.mixin.util.PropertyGroup(properties);
    % end
    % function groups = getPropertyGroups(obj)
    %   properties = struct(...
    %     );
    %   groups = matlab.mixin.util.PropertyGroup(properties);
    % end



  end    % methods(Access=public)
  methods



    % NOTE: Does not work before adding first SDCD.
    function nWholeRowIsNan = get.nWholeRowIsNan(obj)
      nWholeRowIsNan = 0;
      SdcdCa         = obj.Dict.values;

      for i = 1:numel(SdcdCa)
        Sdcd           = SdcdCa{i};
        nWholeRowIsNan = nWholeRowIsNan + sum(Sdcd.bWholeRowIsNan);
      end
    end



  end



end
