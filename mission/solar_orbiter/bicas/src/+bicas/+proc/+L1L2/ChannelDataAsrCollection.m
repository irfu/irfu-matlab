%
% NOTE: Handle class.
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
% NOTE: The constructor does not initialize the entire object completely because
% the constructor call would be too large and awkward then.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ChannelDataAsrCollection < handle
  % PROPOSAL: Better name
  %   ~ASR (as opposed to non-ASR)
  %   ~9x ASR (as opposed to 5x BLTS)
  %   ~9x ASR (as opposed to just one channel)
  %   ~SDID (ASR SDID)
  %   ~SingleChannelData
  %     Which may be used for BLTS(?)
  %   ~Dictionary
  %     CON: Implies arbitrary collection (keys).
  %   ~Collection
  %   ~handle
  %   --
  %   NOTE: SCHD = SingleChannelData
  %   ACHD = AsrChannelData
  %   CDAC = ChannelDataAsrCollection -- IMPLEMENTED
  %   ACHDC, ACDC = AsrChannelDataCollection
  %     PRO: There is need to specify that a channel is an ASR, rather than a
  %          BLTS.
  %     CON: ACHDC is a long abbreviation.
  %     CON: "ACDC" is almost an abbreviation that someone might use though it
  %          has been found in the source code.
  %
  % PROPOSAL: Better tests.
  %
  % PROPOSAL: Implement custom print version of class.  (=??!!!)
  %   PRO: Useful for debugging.
  %
  % PROPOSAL: Constructor must set values for all keys.
  %   PRO: Class is always initialized.
  %   PROPOSAL: Cell array, 9x2, (iSsid, 1)=SSID, (iSsid, 2)=SCHD
  %     CON: Memory problems?
  %   PRO: More natural to assert similarities between constituent SCHDs.
  %
  % PROPOSAL: Make more generic.
  %   CON: Current implementation derives nWholeRowIsNan, which is not generic.
  %     CON-PROPOSAL: Convert to separate function (not method).
  %   PROPOSAL: Store values which are always emulate the size/shape of
  %             1-D column arrays. Values have the same MATLAB class.
  %     Ex: SCHD, FPA, normal arrays.
  %     PRO: Can be used for submitting only VSIB arrays to
  %          bicas.proc.L1L2.qual.get_quality_ZVs_saturation() (without samples).
  %   PROPOSAL: Store objects which always have the same
  %             number of rows (arbitrary size in other dimensions).
  %     PRO: Can store all zVariables, both regular arrays and FPAs.
  %       NOTE: Can then not constrain keys to only be ASR SDIDs.
  %             Need string keys to identify zVariables (hardcoded).
  %     PRO?: Replace bicas.utils.SameRowsMap.
  %     PRO: Can use as argument for
  %          bicas.proc.L1L2.qual.get_QRCBs_channel_saturationn().
  %     PROPOSAL: The object/class itself emulates a 1-D array.
  %       PRO: Class can be used recursively.
  %       PRO: Can extract object for sets of row indices (also recursively).
  %       PRO: Can set object for sets of row indices.
  %         NOTE: Cf. bicas.utils.SameRowsMap.set_rows()
  %
  % PROPOSAL: Assert that all SCHDs have the same size (columns) and data types.



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
  properties(GetAccess=public, SetAccess=private)
    nRecords
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Initialize empty object, without any channels.
    function obj = ChannelDataAsrCollection(nRecords)
      assert(isscalar(nRecords) & isnumeric(nRecords) & nRecords >= 0)

      obj.Dict     = configureDictionary('uint8', 'cell');
      obj.nRecords = nRecords;
    end



    % Set channel data. Either (1) add new channel, or (2) overwrite
    % pre-existing channel.
    function obj = set_channel(obj, asrSdid, Schd)
      % IMPLEMENTATION NOTE: Can not assert that the caller does not overwrite
      % the CHSD for a specified ASR SDID since overwriting is required by
      % bicas.proc.L1L2.demuxer.reconstruct_ASR_samples_NEW()
      % /reconstruct_missing_data_helper() (inner function).
      assert(isscalar(asrSdid))
      assert(ismember(asrSdid, bicas.proc.L1L2.const.C.SDID_ASR_AR))
      assert(isa(Schd, 'bicas.proc.L1L2.SingleChannelData'))
      assert(size(Schd, 1) == obj.nRecords)

      obj.Dict(asrSdid) = {Schd};
    end



    function Schd = get_channel(obj, asrSdid)
      assert(isscalar(asrSdid))
      ca   = obj.Dict(asrSdid);
      Schd = ca{1};
    end



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
