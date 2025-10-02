%
% Class that represents how to convert particular combinations of QRCID and DSI
% into modifications of quality ZVs, data blanking etc. for multiple DSIs.
% Stores up to one QRCS per QRCID+DSI.
%
% NOTE: Is handle class.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingsMap < handle & matlab.mixin.Copyable
  % PROPOSAL: Rename QrcsMap.
  %   PRO: More analogous to QrcbMap.
  %
  % PROPOSAL: Require that all QRCSs are of a certain subclass. Add MC to
  %           constructor.
  %
  % PROBLEM: The same QRCID may apply to both L1/L1R-->L2 and L2-->L3 processing,
  %          and may be read from NSO table for both.
  %   Ex: Floating ANT3 after sweeps labelled in L2, not used for L3 EFIELD.
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
  %       Ex: Floating ANT3 after sweeps.
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
  %
  % ============================================================================
  % NOTE: Important QRC functions for QRC settings:
  %   bicas.proc.qrc.NSO_table_to_QRCBM(requestedQrcidAr, NsoTable, tt2000Ar, L))
  %     NOTE: No QRCSs needed, except QRCIDs for the type of processing.
  %   bicas.proc.L2L3.qrc.L2QBM_to_QRCBs(l2qbmAr, lxqbmName, Qrcsm)
  %     NOTE: Should ideally only use the subset of QRCSs for which quality bits
  %           are read.
  %       NOTE: Uses QRCSs for processing which created the quality bits
  %             (L1/L1R-->L2), but sets QRCBs for QRCIDs used in the processing
  %             that will use the information (L2-->L3).
  %   bicas.proc.qrc.QRCB_arrays_to_quality_ZVs(Qrcbm, Qrcsm, lxqbmName))
  %   bicas.proc.L1L2.qrc.set_5xBLTS_voltage_samples_FV(voltageAr, ssidAr, Qrcbm, Qrcsm)
  %   bicas.proc.L1L2.qrc.set_current_samples_FV(currentAr, Qrcbm, Qrcsm)
  %   bicas.proc.L2L3.qrc.L2QBM_to_channel_saturation_QRCBs(l2qbmAr, saturationQualitySchemeId)
  %   bicas.proc.L2L3.qrc.set_FPA_samples_FP(Fpa, Qrcbm, Qrcsm, qrcsFieldName)
  %
  % NOTE: QRC functions which do not need QRCSs but do need subsets of QRCs
  %       (subset relative to a type of processing):
  % bicas.proc.L1L2.qrc.VSIBs_to_saturation_QRCBs(tt2000Ar, saturationQualitySchemeId, VsibZvm, isSwf, vstbFractionThreshold, cwfSlidingWindowLengthSec)
  %   Iterates over CHANNEL_SATURATION QRCIDs.
  %   Calls bicas.proc.L1L2.qrc.channel_saturation_to_global_saturation_QRCBs(ChannelSaturationQrcbm, nRecords).
  % bicas.proc.L1L2.qrc.channel_saturation_to_global_saturation_QRCBs(ChannelSaturationQrcbm)
  %   Iterates over CHANNEL_SATURATION QRCIDs.
  % bicas.proc.L2L3.qrc.L2QBM_to_channel_saturation_QRCBs(l2qbmAr, saturationQualitySchemeId)
  %   Wrapper around bicas.proc.L2L3.qrc.L2QBM_to_QRCBs().
  %   NOTE: Needs set of QRCSs for which quality bits to read
  %         (CHANNEL_SATURATION QRCIDs).
  % ============================================================================



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Dependent)
    qrcidAr
  end
  properties(SetAccess=private, GetAccess=private)
    QrcsDict
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods



    % Return list of QRCIDs.
    function qrcidAr = get.qrcidAr(obj)
      % IMPLEMENTATION NOTE: Sort to gain deterministic result which is good for
      % testing.
      qrcidAr = sort(obj.QrcsDict.keys);
      assert(isstring(qrcidAr) & iscolumn(qrcidAr))
    end



  end
  methods(Access=public)



    function obj = QrcSettingsMap()
      % IMPLEMENTATION NOTE: configureDictionary("cell",
      % "bicas.proc.QrcSetting") does not work since dictionary apparently does
      % not permit specifying abstract classes. Therefore using cell values.
      % MATLAB error message: "Unable to configure dictionary with abstract
      % type 'bicas.proc.QrcSetting'."

      obj.QrcsDict = configureDictionary("string", "cell");
    end



    % Add new QRCS. Must not reuse QRCID (assertion).
    function add(obj, qrcid, Qrcs)
      assert(isa(Qrcs, "bicas.proc.QrcSetting") & isscalar(Qrcs))
      obj.assert_scalar_QRCID(qrcid);

      assert(~obj.QrcsDict.isKey(qrcid))

      obj.QrcsDict(qrcid) = {Qrcs};   % NOTE: Cell in cell.
    end



    function remove(obj, qrcid)
      assert(isstring(qrcid) & isscalar(qrcid))
      assert(obj.QrcsDict.isKey(qrcid))

      obj.QrcsDict = obj.QrcsDict.remove(qrcid);
    end



    function remove_many(obj, qrcidAr)
      assert(iscolumn(qrcidAr))
      for qrcid = qrcidAr'
        obj.remove(qrcid)
      end
    end



    % RETURN VALUE
    % ============
    % QRCS for specified key. Key must exist.
    function Qrcs = get(obj, qrcid)
      obj.assert_scalar_QRCID(qrcid);
      assert(obj.QrcsDict.isKey(qrcid), ...
        "There is no key qrcid=""%s"".", qrcid)

      QrcsCa = obj.QrcsDict(qrcid);
      Qrcs   = QrcsCa{1};
    end



    % NOTE: Asserts against QRCID collisions.
    function add_QRCSM(obj, Qrcsm)
      assert(isa(Qrcsm, "bicas.proc.QrcSettingsMap"))

      % IMPLEMENTATION NOTE: dictionary is not a handle object.
      TempDict = obj.QrcsDict;

      % Merge dictionaries in temporary copy (so can abort in case of failure).
      TempDict(Qrcsm.QrcsDict.keys) = Qrcsm.QrcsDict.values;

      % ASSERT: No key collisions.
      assert(TempDict.numEntries == ...
        obj.QrcsDict.numEntries + Qrcsm.QrcsDict.numEntries, ...
        "QRCSMs have overlapping QRCIDs.")

      obj.QrcsDict = TempDict;
    end



  end    % methods(Access=public)



  %############################
  %############################
  % PROTECTED INSTANCE METHODS
  %############################
  %############################
  methods(Access = protected)



    % Support deep copies with copy() by overriding matlab.mixin.copyable
    % method.
    function QrcbmCopy = copyElement(obj)
      QrcbmCopy = bicas.proc.QrcSettingsMap();

      for qrcid = obj.qrcidAr'
        QrcbmCopy.add(qrcid, obj.get(qrcid))
      end
    end



  end



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % NOTE: Assert scalar QRCID, in particular to prevent using an array of
    % keys on the internal dictionary.
    function assert_scalar_QRCID(qrcid)
      assert(isstring(qrcid) & isscalar(qrcid))
    end



  end    % methods(Static, Access=private)



end
