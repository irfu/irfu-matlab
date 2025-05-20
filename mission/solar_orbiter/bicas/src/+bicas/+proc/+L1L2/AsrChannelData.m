%
% NOTE: Handle class
%
% Class which effectively wraps a dictionary
% (ASR) SDID-->bicas.proc.L1L2.SingleChannelData.
%
% This useful since
% (1) (MATLAB) dictionary values can not be non-scalar
%     (bicas.proc.L1L2.SingleChannelData is a column vector). Therefore,
%     the implementation can work around this.
% (2) it can sum up the number of SCHD fill positions.
%
% NOTE: The constructor does not initialize the object completely because
% the constructor call would be too large and awkward then.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef AsrChannelData < handle
  % PROPOSAL: Better name
  %   ~ASR (as opposed to non-ASR)
  %   ~9x ASR (as opposed to 5x BLTS)
  %   ~9x ASR (as opposed to just one channel)
  %   ~SingleChannelData
  %     Which may be used for BLTS(?)
  %   ~Dictionary
  %     CON: Implies arbitrary collection (keys).
  %   ~Collection
  %   ~handle
  %   --
  %   ASCHD = AsrSingleChannelData
  %     CON: Does not imply collection of ALL ASR channels.
  %       PRO: "Single" in the name.
  %   ACHD = AsrChannelData -- IMPLEMENTED
  %   ACHDC, ACDC = AsrChannelDataCollection
  %     PRO: There is need to specify that a channel is an ASR, rather than a
  %          BLTS.
  %     CON: ACHDC is a long abbreviation.
  %     CON: "ACDC" is almost an abbreviation that someone might use though it
  %          has been found in the source code.
  %
  % PROPOSAL: Better tests.
  %
  % PROPOSAL: Constructor pre-allocates SCHDs.
  % PROPOSAL: Implement custom print version of class.
  %   PRO: Useful for debugging.
  %
  % PROPOSAL: Constructor must set values for all keys.
  %   PRO: Class is always initialized.
  %   PROPOSAL: Cell array, 9x2, (iSsid, 1)=SSID, (iSsid, 2)=SCHD
  %     CON: Memory problems?
  %   PRO: More natural to assert similarities between constituent SCHDs.
  %
  % PROPOSAL: Assert that all SCHDs have the same size and data types.



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Dependent)
    % Total number of bWholeRowIsNan values in all the SCHD objects stored
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



    % Initialize empty object, without any channels.
    function obj = AsrChannelData()
      obj.Dict = dictionary;
    end



    % Set channel data.
    function obj = set_channel(obj, asrSdid, Schd)
      assert(isscalar(asrSdid))
      assert(ismember(asrSdid, bicas.proc.L1L2.const.C.SDID_ASR_AR))
      assert(isa(Schd, 'bicas.proc.L1L2.SingleChannelData'))

      obj.Dict(asrSdid) = {Schd};
    end



    function Schd = get_channel(obj, asrSdid)
      assert(isscalar(asrSdid))
      ca   = obj.Dict(asrSdid);
      Schd = ca{1};
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



    % NOTE: Does not work before adding first SCHD.
    function nWholeRowIsNan = get.nWholeRowIsNan(obj)
      nWholeRowIsNan = 0;
      SchdCa         = obj.Dict.values;

      for i = 1:numel(SchdCa)
        Schd           = SchdCa{i};
        nWholeRowIsNan = nWholeRowIsNan + sum(Schd.bWholeRowIsNan);
      end
    end



  end



end
