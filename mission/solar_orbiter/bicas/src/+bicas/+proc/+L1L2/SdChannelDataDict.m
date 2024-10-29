%
% Class which effectively wraps a dictionary
% SDID-->bicas.proc.L1L2.SdChannelData.
% This useful since
% (1) (MATLAB) dictionary values can not be non-scalar
%     (bicas.proc.L1L2.SdChannelData is a column vector). Therefore,
%     the implementation must work around this (it uses 1x1 cell arrays).
% (2) it can sum up the number of SDCD fill positions.
%
% NOTE: The constructor does not initialize the object completely (because
% the constructor call would be too awkward).
%
% NOTE: Is NOT a handle class. Should maybe become, if performance becomes a
% problem?
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SdChannelDataDict
  % PROPOSAL: Automatic test code.
  %
  % PROPOSAL: Constructor pre-allocates SDCDs.
  % PROPOSAL: Implement custom print version of class.
  %   PRO: Useful for debugging.



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant, Access=private)
    KEYS_AR = bicas.proc.L1L2.const.C.SDID_ASR_AR
  end
  properties(Dependent)
    % Total number of fill positions in the underlying
    nFp
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



    function obj = set(obj, asrSdid, Sdcd)
      assert(isscalar(asrSdid))
      assert(isa(Sdcd, 'bicas.proc.L1L2.SdChannelData'))
      assert(ismember(asrSdid, obj.KEYS_AR))
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



    function nFp = get.nFp(obj)
      nFp    = 0;
      SdcdCa = obj.Dict.values;
      for i = 1:numel(SdcdCa)
        Sdcd = SdcdCa{i};
        nFp  = nFp + sum(Sdcd.bFp);
      end
    end



  end



end
