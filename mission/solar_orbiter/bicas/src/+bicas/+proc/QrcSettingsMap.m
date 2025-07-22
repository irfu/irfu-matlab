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
  % PROBLEM: Misspelled DSIs leads get() to assume that it is a DSI for which
  %          there is no QRCS.
  %   PROPOSAL: QrcSettingsMap constructor accepts list of all possible DSIs and asserts
  %             using it.
  %
  %
  %
  % PROBLEM: The same QRCID may apply to both L1R-->L2 and L2-->L3 processing,
  %          may be read from NSO table. -- ~IMPLEMENTED
  %   Ex: Floating V3 after sweeps labelled in L2, not used for L3 EFIELD.
  %       https://github.com/irfu/irfu-matlab/issues/156
  %   CON: QRCs can apply to L2-->L3 but not L1R-->L2, which makes them very
  %        artificial being used in L1R-->L2 code.
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
  %   --
  %   PROPOSAL: Have QRCSs represent *all* types of processing. -- IMPLEMENTED
  %     PROPOSAL: QRCS (QrcSettings) owns separate "sub-QRCS" objects/classes
  %               for different kinds of processing.
  %       PRO: Can reference "entire set of QRC settings" for type of
  %            processing.
  %       PROPOSAL: Singleton "database" which is map from QRCID+DSI-->~QRCS.
  %         NOTE: Assumes that there is no need for information associated with
  %               QRCID only.
  %       PROPOSAL: Have one sub-QRCS per DSI.
  %         CON?: How handle different types of processing producing the same
  %               DSI?
  %           CON: Should only happen for unofficial processing.
  %         CON: Must be initialized with many default values.
  %           CON-PROPOSAL: Absent sub-QRCS implies defaults.
  %             PRO: Can handle multiple sub-QRCS classes.
  %         TODO-DEC: Class names, abbreviations?
  %           QRCS=QrcSetting: All settings for QRCID
  %           QRCS=QrcSetting: All settings for one QRCID+DSI.
  %       PROPOSAL: bicas.proc.qual.QRCB_arrays_to_quality_ZVs() should have
  %                 argument for (one) DSI. Only execute the needed DSI.
  %
  %
  %
  % PROPOSAL: Treat input L2QBM as another generic (though not completely)
  %           input source for QRCBs, next to NSO table and autodetection,
  %           using QRCS L2QBM bit.
  %   PROBLEM: How distinguish between QRCs which can be determines this way
  %            and not?
  %     PROPOSAL: There is exactly one L2QBM quality bit.
  %     PROPOSAL: Enumerate QRCIDs (hardcoded).
  %     PROPOSAL: Assert that QRCS defines exactly one (unique?) L2QBM quality
  %               bit.
  %     PROPOSAL: Dedicated flag in QRCS.
  %     CON-PROPOSAL: Only implement for QRCs for which it is relevant.
  %                   Having general system is unnecessary, too broad.
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
  % PROPOSAL: QRCSs should not be function of DSI, but some
  %           ~"processing type ID" which can group output DSIs arbitrarily.
  %   PROPOSAL: Replace QRCS-->QRCGS=QrcGroupSettings.
  %
  %
  %
  % TODO: Remove QrcSettings. Eliminate QRCS abbreviation. -- DONE
  % TODO: Rename QrcSettingsDb --> QrcsMap=QrcSettingsMap(?) -- DONE
  % TODO: Convert QrcDsiSettings --> ~QRCS=QrcSettings? -- DONE
  % TODO: DSI arguments--> ~ProcessingType ID?
  %
  % TODO-DEC: Require legal DSI/PTID when initializing?
  % TODO-DEC: Require existing QRCID+DSI/PTID combination when calling?
  %   PROPOSAL: Method for "locking" object. After that point, .get() always
  %             verifies against the separate lists of QRCIDs and DSIs found in
  %             the object.
  %     CON: May work for QRCIDs, but not necessarily for DSIs. Not all DSIs
  %          might have any QRCS.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private, GetAccess=private)
    QrcsDict
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = QrcSettingsMap()
      % IMPLEMENTATION NOTE: configureDictionary("cell",
      % "bicas.proc.QrcSetting") does not work since dictionary apparently
      % does not permit abstract classes. Therefore using cell values.
      % MATLAB error message: "Unable to configure dictionary with abstract type
      % 'bicas.proc.QrcSetting'."
      obj.QrcsDict = configureDictionary("cell", "cell");
    end



    % Add new QRCS. Can not overwrite (assertion).
    function add(obj, qrcid, dsiCa, Qrcs)
      assert(isa(Qrcs, "bicas.proc.QrcSetting") & isscalar(Qrcs))
      assert(iscell(dsiCa) & iscolumn(dsiCa))

      for dsiCa0 = dsiCa'
        dsi = dsiCa0{1};
        key = obj.assert_get_dictionary_key(qrcid, dsi);

        assert(~obj.QrcsDict.isKey(key))

        obj.QrcsDict(key) = {{Qrcs}};   % NOTE: Cell in cell.
      end
    end



    % RETURN VALUE
    % ============
    % QRCS for specified key, if key exists.
    % Empty, if key does not exist.
    function Qrcs = get(obj, qrcid, dsi)
      key = obj.assert_get_dictionary_key(qrcid, dsi);

      if obj.QrcsDict.isKey(key)
        QrcidCaCa = obj.QrcsDict(key);
        Qrcs     = QrcidCaCa{1}{1};
      else
        Qrcs = [];
      end
    end



    % Return list of unique QRCIDs.
    function qrcidAr = get_QRCIDs(obj)
      qrcidAr = string.empty(0, 1);

      for keyCaCa = obj.QrcsDict.keys'
        qrcid = keyCaCa{1}{1};

        qrcidAr(end+1, 1) = qrcid;
      end

      % NOTE: unique() should convert column-->column.
      qrcidAr = unique(qrcidAr);
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

      obj.QrcsDict = TempDict;
    end



  end    % methods(Access=public)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)
  end    % methods(Access=private)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)
  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function mapKey = assert_get_dictionary_key(qrcid, dsi)
      assert(isstring(qrcid))
      assert(ischar(dsi))

      % NOTE: Must use cell array within cell array. Dictionary will vectorize
      % the outer cell array.
      % IMPLEMENTATION NOTE: containers.Map() will not allow using
      % non-vectorized cell arrays at all.
      mapKey = {{qrcid; dsi}};
    end



  end    % methods(Static, Access=private)



end
