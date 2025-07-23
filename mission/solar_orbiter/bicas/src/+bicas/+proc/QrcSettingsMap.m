%
% Class that represents how to convert particular combinations of QRCID and DSI
% into modifications of quality ZVs, data blanking etc. for multiple DSIs.
% Stores up to one QRCS per QRCID+DSI.
%
% NOTE: Is handle class because of convenience.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingsMap < handle
  % PROPOSAL: Rename QrcsMap.
  %   PRO: More analogous to QrcbMap.
  %
  % PROPOSAL: Change get_QRCIDs() to dependent property qrcidAr.
  %
  %
  % PROBLEM: The same QRCID may apply to both L1/L1R-->L2 and L2-->L3 processing,
  %          and may be read from NSO table for both.
  %   Ex: Floating V3 after sweeps labelled in L2, not used for L3 EFIELD.
  %       https://github.com/irfu/irfu-matlab/issues/156
  %   Ex: ANTx_FAIL
  %   NOTE: QRCs may also ONLY apply to L2-->L3 and NOT L1R-->L2, which makes
  %         them very artificial being used in L1R-->L2 code.
  %     Ex: Bad density.
  %   --
  %   PROPOSAL: Read quality bits set in L2 when processing L2-->L3.
  %     PRO: State in L2 is always in sync with L3 when processing.
  %     CON: Must interpret quality bits and use for setting QRCBs (or similar)
  %          for L2-->L3. System is not designed for that.
  %       CON-PROPOSAL: Design standard functionality for setting QRCBs from
  %                     quality bits (L2QBM), using QRCSs(?).
  %   --
  %   PROPOSAL: Use NSO table to initialize QRCBs for both L1R->L2 and L2->L3
  %             processing (same QRCIDs). ==> Every QRCS contains settings for
  %             both L1R->L2 and L2->L3 (and technically all other forms of
  %             processing, e.g. L2-->L2).
  %     PRO: Information does not need to be unambiguously available in quality
  %          zVariables.
  %       Ex: Floating V3 after sweeps.
  %       Ex: QUALITY_FLAG is lowered, but it is not known for which reason:
  %         Ex: Other effect, BICAS's global cap, ROC's global cap.
  %     PRO: Independent of how information is represented in zVariables.
  %       Ex: Bit locations, bit combinations.
  %         Ex: Full/partial saturation represented through two bits. Must read
  %             both to determine whether full+partial or only partial
  %             saturation is set.
  %     CON: Only works for QRCBs which are set through NSO table, not
  %          autodetected.
  %       Ex: Saturation: NSO table+autodetected (nominal)
  %     PROBLEM: QRCs must store both L2QBM and L3QBM which creates problem for
  %              code that is generic to LxQBM.
  %       PROBLEM: Need to distinguish L2QBM for L1R-->L2 and L2-->L2.
  %     CON-PROPOSAL: Use multiple inheritance, but only when relevant.
  %
  %
  %
  % PROBLEM: How represent that L3 density can raise the QUALITY_FLAG compared
  %          to parent?
  %   NOTE: The only case is when derive L3 (density?) without saturated
  %         channels.
  %   NOTE: Can only derive higher QUALITY_FLAG if can derive what L2(!)
  %         QUALITY_FLAG should have been without saturation?! ==> Must do this
  %         L2 derivation when processing L2-->L3?!!
  %
  %
  %
  % PROBLEM: How hardcode QRCSs for L3 OSR and DSR respectively?
  %   PROPOSAL: Implicit that L3 OSR QRCSs apply also to DSRs?
  % PROBLEM: How hardcode blanking of input data to L3?
  %          L2 CWF-->L3 EFIELD+SCPOT; L3 SCPOT-->L3 DENSITY.
  %   PROPOSAL: See L3 DENSITY as function of L3 SCPOT (CDF to CDF).
  %             ==> Implicit that QRCS for L3 SCPOT affects L3 DENSITY.
  %     PROBLEM: Still one variable for both L3 EFIELD+SCPOT.
  %
  %
  % PROPOSAL: Convert class to only maintaining map QRCID-->QRCS, for QRCSs for
  %           ONLY the same type of processing. Use separate global constants for
  %           storing different subsets of QRCSs, e.g. for different types of
  %           processing. If needs different QRC settings for e.g. different
  %           DSIs for the same type of processing, then that should be selected
  %           outside of the qual functions which use QRC settings.
  %   CON: Destroys any link between the same QRCID for different types of
  %        processing.
  %     CON: Code execution assures some correspondence (should): identical
  %          QRCIDs.
  %     CON: Code which assigns maps should co-locate related QRCSs.
  %
  % PROPOSAL: Require that all QRCSs are of a certain subclass. Add MC to
  %           constructor.
  %
  %
  %
  % PROPOSAL: qual functions should only need to know about mapping
  %   QRCID-->action (e.g. QRCS). They should not need to know how the QRCSs
  %   were selected, e.g. due to PTID or DSI.
  % ============================================================================
  % NOTE: Important qual functions for QRC settings:
  %   bicas.proc.qual.NSO_table_to_QRCB_map(requestedQrcidAr, NsoTable, tt2000Ar, L))
  %     NOTE: No QRCSs needed, except QRCIDs for the type of processing.
  %   bicas.proc.qual.LxQBM_to_QRCB_maps(l2qbmAr, lxqbmName, QrcsMap)
  %     EXPERIMENTAL. UNUSED.
  %     NOTE: Should ideally only use the subset of QRCSs for which quality bits
  %           are read.
  %       NOTE: Uses QRCSs for processing which created the quality bits
  %             (L1/L1R-->L2), but sets QRCBs for QRCIDs used in the processing
  %             that will use the information (L2-->L3).
  %   bicas.proc.qual.QRCB_arrays_to_quality_ZVs(QrcbMap, Qrcsm, ptid, lxqbmName))
  %   bicas.proc.L1L2.qual.set_5xBLTS_voltage_samples_FV(voltageAr, ssidAr, QrcbMap, Qrcsm, ptid)
  %   bicas.proc.L1L2.qual.set_current_samples_FV(currentAr, QrcbMap, Qrcsm, ptid)
  %   bicas.proc.L2L3.get_saturation_QRCBs_from_L2QBM(l2qbmAr, saturationQualitySchemeId)
  %     EXPERIMENTAL. UNUSED.
  %   bicas.proc.L2L3.set_VDC_EDC_samples_FV(VDC_Fpa, EDC_Fpa, QrcbMap, QrcsMap)
  %     EXPERIMENTAL. UNUSED.
  %
  % NOTE: qual functions which do not need QRCSs but do need subsets of QRCs (subset relative to a type of
  %       processing):
  % bicas.proc.L1L2.get_saturation_QRCBs(tt2000Ar, saturationQualitySchemeId, VsibZvm, isSwf, vstbFractionThreshold, cwfSlidingWindowLengthSec)
  %   Iterates over CHANNEL_SATURATION QRCIDs.
  %   Calls bicas.proc.L1L2.qual.channel_saturation_to_global_saturation_QRCBs(ChannelSaturationQrcbMap, nRecords).
  % bicas.proc.L1L2.qual.channel_saturation_to_global_saturation_QRCBs()
  %   Iterates over CHANNEL_SATURATION QRCIDs.
  % bicas.proc.L2L3.get_saturation_QRCBs_from_L2QBM(l2qbmAr, saturationQualitySchemeId)
  %   EXPERIMENTAL. UNUSED.
  %   Wrapper around bicas.proc.qual.LxQBM_to_QRCB_maps().
  %   NOTE: Needs set QRCSs for which quality bits to read (CHANNEL_SATURATION QRCIDs).
  % ============================================================================



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private, GetAccess=private)
    QrcsDict
    % legalPtidAr
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % legalPtidAr
    %       Column array of string constants. These are the only allowed PTIDs
    %       that can be used with this instance.
    %       RATIONALE: This is a safeguard against typos.
    %
    function obj = QrcSettingsMap()
      % IMPLEMENTATION NOTE: configureDictionary("cell",
      % "bicas.proc.QrcSetting") does not work since dictionary apparently
      % does not permit abstract classes. Therefore using cell values.
      % MATLAB error message: "Unable to configure dictionary with abstract type
      % 'bicas.proc.QrcSetting'."

      % assert(isstring(legalPtidAr) & iscolumn(legalPtidAr))

      obj.QrcsDict    = configureDictionary("string", "cell");
      % obj.legalPtidAr = legalPtidAr;
    end



    % Add new QRCS. Can not overwrite (assertion).
    function add(obj, qrcid, Qrcs)
      assert(isa(Qrcs, "bicas.proc.QrcSetting") & isscalar(Qrcs))

      key = obj.assert_get_dictionary_key(qrcid);

      assert(~obj.QrcsDict.isKey(key))

      obj.QrcsDict(key) = {Qrcs};   % NOTE: Cell in cell.
    end



    % RETURN VALUE
    % ============
    % QRCS for specified key. Key must exist.
    function Qrcs = get(obj, qrcid)
      key = obj.assert_get_dictionary_key(qrcid);

      assert(obj.QrcsDict.isKey(key), ...
        "There is no key qrcid=""%s"".", qrcid)

      QrcsCa = obj.QrcsDict(qrcid);
      Qrcs   = QrcsCa{1};

      % assert(obj.QrcsDict.isKey(key), ...
      %   "There is no key for qrcid=""%s"" and ptid=""%s"".", qrcid, ptid)
      % if obj.QrcsDict.isKey(key)
      %   QrcidCaCa = obj.QrcsDict(key);
      %   Qrcs     = QrcidCaCa{1}{1};
      % else
      %   Qrcs = [];
      % end
    end



    % Return list of unique QRCIDs.
    function qrcidAr = get_QRCIDs(obj)
      % qrcidAr = string.empty(0, 1);
      %
      % for keyCaCa = obj.QrcsDict.keys'
      %   qrcid = keyCaCa{1}{1};
      %
      %   qrcidAr(end+1, 1) = qrcid;
      % end
      %
      % % NOTE: unique() should convert column-->column.
      % qrcidAr = unique(qrcidAr);

      % IMPLEMENTATION NOTE: Sort to gain deterministic result which is good for
      % testing.
      qrcidAr = sort(obj.QrcsDict.keys);
      assert(isstring(qrcidAr) & iscolumn(qrcidAr))
    end



    function merge(obj, Qrcsm)
      assert(isa(Qrcsm, "bicas.proc.QrcSettingsMap"))

      % IMPLEMENTATION NOTE: dictionary is not a handle object.
      TempDict = obj.QrcsDict;

      % Merge dictionaries in temporary copy (so can abort in case of failure).
      TempDict(Qrcsm.QrcsDict.keys) = Qrcsm.QrcsDict.values;

      % ASSERT: No key collisions.
      assert(TempDict.numEntries == ...
        obj.QrcsDict.numEntries + Qrcsm.QrcsDict.numEntries)

      obj.QrcsDict    = TempDict;
      % obj.legalPtidAr = unique([obj.legalPtidAr; Qrcsm.legalPtidAr]);
    end



    % RETURN VALUE
    % ============
    % containers.Map QRCID-->QRCS for the specified PTID.
    %
    % function QrcsMap = get_QRCID_QRCS_map(obj, ptid)
    %   assert(ismember(ptid, obj.legalPtidAr))
    %
    %   QrcsMap = containers.Map();
    %
    %   for key = obj.QrcsDict.keys'
    %     currentQrcid = key{1}{1};
    %     currentPtid  = key{1}{2};
    %     QrcsCaCa    = obj.QrcsDict(key);
    %     Qrcs        = QrcsCaCa{1}{1};
    %
    %     if ptid == currentPtid
    %       QrcsMap(currentQrcid) = Qrcs;
    %     end
    %   end
    % end



  end    % methods(Access=public)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function mapKey = assert_get_dictionary_key(obj, qrcid)
      % PROPOSAL: Abolish. Exists for historical reasons only.

      assert(isstring(qrcid) & isscalar(qrcid))
      % assert(isstring(ptid) & isscalar(ptid))
      % assert(ismember(ptid, obj.legalPtidAr))

      % NOTE: Must use cell array within cell array. Dictionary will vectorize
      % the outer cell array.
      % IMPLEMENTATION NOTE: containers.Map() will not allow using
      % non-vectorized cell arrays at all.
      % mapKey = {{qrcid; ptid}};
      mapKey = qrcid;
    end



  end    % methods(Access=private)



end
