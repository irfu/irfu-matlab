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
