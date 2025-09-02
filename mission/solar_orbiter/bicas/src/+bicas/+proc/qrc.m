%
% Collection of reusable, generic code relating to QRCs.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qrc



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Obtain QRCBs from NSO table for specified timestamps.
    %
    %
    % ARGUMENTS
    % =========
    % requestedQrcidAr
    %       Column array of QRCIDs for which the return value shall contain
    %       data. This is both for being able to
    %       (1) filter the content of the NSO table (only return data for
    %           selected QRCIDs), and
    %       (2) return default values (QRCB=false) for QRCIDs for which there
    %           are no events which overlap with the specified time stamps (or
    %           no such events at all).
    %
    %
    % RETURN VALUE
    % ============
    % Qrcbm
    %       Contains keys for all QRCIDs specified in requestedQrcidAr, not just
    %       those present in NsoTable.
    %
    function Qrcbm = NSO_table_to_QRCBM(...
        requestedQrcidAr, NsoTable, tt2000Ar, L)

      % IMPLEMENTATION NOTE: Without requestedQrcidAr, the function can not
      % create a return value map that contains keys for all QRCIDs, in case
      % the NsoTable does not contain all QRCIDs.

      % Local variable naming conventions:
      % ----------------------------------
      % GE = Global Event = NSO event in global NSO event table.
      % LE = Local Event  = NSO event that overlaps with the specified
      %                     timestamps.

      assert(isstring(requestedQrcidAr) & iscolumn(requestedQrcidAr))

      NsoEventMatchAr = NsoTable.get_NSO_event_matches(tt2000Ar);

      nLe = numel(NsoEventMatchAr);
      nGe = NsoTable.nEvents;
      L.logf('info', ...
        ['Searched non-standard operations (NSO) table.', ...
        ' Found %i relevant NSO events out of a total of %i NSO events.'], ...
        nLe, nGe);

      %----------------------------------------------------------------------
      % Initialize "empty" Qrcbm (elements=false) for all requested QRCIDs
      %----------------------------------------------------------------------
      % IMPLEMENTATION NOTE: valueType=logical implies scalar (sic!) and can
      %                      therefore not be used.
      Qrcbm = bicas.proc.QrcbMap(numel(tt2000Ar));
      for i = 1:numel(requestedQrcidAr)
        Qrcbm.add(requestedQrcidAr(i), false(size(tt2000Ar)));
      end

      %-----------------------------------------------------------------------
      % Iterate over NSO events and set QRCBs for resp. QRCIDs and timestamps
      %-----------------------------------------------------------------------
      for kLe = 1:nLe
        eventQrcid = NsoEventMatchAr(kLe).qrcid;

        if ismember(eventQrcid, requestedQrcidAr)
          % Index into GLOBAL NSO events table.
          iGe    = NsoEventMatchAr(kLe).iNsoEvent;
          % Logical indices into tt2000Ar.
          qrbcAr = NsoEventMatchAr(kLe).qrcbAr;

          %=====================================================================
          % Log relevant NSO events by referring to the GLOBAL NSO events table
          %=====================================================================
          L.logf('info', '    %s -- %s %s', ...
            bicas.utils.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGe), 9), ...
            bicas.utils.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array( iGe), 9), ...
            eventQrcid);

          %====================================
          % Update corresponding QRCB elements
          %====================================
          Qrcbm.set(eventQrcid, Qrcbm.get(eventQrcid) | qrbcAr);
        end
      end    % for
    end



    % Given QRCB arrays, translate them into
    % (1) ZV QUALITY_FLAG, and
    % (2) ZV L*_QUALITY_BITMASK.
    %
    % NOTE: Does not work with FPAs, since the internal algorithm can not
    % produce unknown values. The caller is supposed to decide how to interpret
    % unknown values (QRCBs) before calling the function.
    %
    %
    % ARGUMENTS
    % =========
    % Qrcbm
    % Qrcsm
    %       Must contain the same keys as Qrcbm.
    % lxqbmName
    %       String constant/ZV name which represents the QRCS field that should
    %       be referenced for specifying the LxQBM value.
    %
    %
    % RETURN VALUES
    % =============
    % QUALITY_FLAG
    %       QUALITY_FLAG max value wrt. to QRCs handled in this function.
    % lqzbm
    %       L*_QUALITY_BITMASK value wrt. to QRCs handled in this function.
    %       Refers to L2_QUALITY_BITMASK or L3_QUALITY_BITMASK depending on
    %       context.
    %
    function [qfl, lxqbm] = QRCB_arrays_to_quality_ZVs(...
        Qrcbm, Qrcsm, lxqbmName)
      % PROPOSAL: Split into separate functions for QUALITY_FLAG and LxQBM.
      %   PRO: EFIELD and SCPOT do not have L3QBM.
      %   PRO: Simpler-ish testing
      %   CON: More code. Functions will almost be the same.
      % PROPOSAL: Somehow abolish the "lxqbmName" argument.
      %   PRO: Is in principle unnecessary since the QRCS class should imply
      %        it.
      %   PROPOSAL: All (applicable) QRCS classes use the same field name.
      %   PROPOSAL: This function determines the field name from the QRCS class
      %             (using hardcoded table).
      %     CON: Less general.
      %     NOTE: Already uses hardcoded table.

      % Dictionary for translations from "string constant"/ZV name representing
      % a field value, to the actual field value.
      DICT_LXQBM_NAME_TO_FIELD_NAME = dictionary(...
        ["L2_QUALITY_BITMASK", "L3_QUALITY_BITMASK"], ...
        ["l2qbm",              "l3qbm"]);

      assert(isa(Qrcbm, "bicas.proc.QrcbMap"))
      assert(isa(Qrcsm, "bicas.proc.QrcSettingsMap"))
      assert(isequal(Qrcbm.qrcidAr, Qrcsm.qrcidAr))
      assert(isstring(lxqbmName))

      nRec           = Qrcbm.nRecords;
      lxqbmFieldName = DICT_LXQBM_NAME_TO_FIELD_NAME(lxqbmName);

      % Create "empty" quality variable arrays, with max possible quality
      % (QUALITY_FLAG max, Lx_QUALITY_BITMASK=0), which can then later be
      % "lowered" if necessary.
      qfl   = ones( nRec, 1, 'uint8' ) * bicas.const.qrc.QFL_MAX;
      lxqbm = zeros(nRec, 1, 'uint16');

      %=================================
      % Iterate over QRCIDs in argument
      %=================================
      for qrcid = Qrcbm.qrcidAr'
        qrcbAr = Qrcbm.get(qrcid);
        Qrcs   = Qrcsm.get(qrcid);

        %------------------
        % Set QUALITY_FLAG
        %------------------
        % IMPLEMENTATION NOTE: Only adjusts relevant indices since the
        % operation is more natural (simpler, shorter) that way. min() should
        % only be applied to the indices where QRCB=true.
        qfl(qrcbAr) = min(qfl(qrcbAr), Qrcs.qfl);

        %--------------
        % Update LXQBM
        %--------------
        lxqbm = bitor(...
          lxqbm, ...
          Qrcs.(lxqbmFieldName) * uint16(qrcbAr));
      end
    end



    % Given QRCB arrays, translate them into GA "CAVEATS".
    %
    % NOTE: This function does not (and should not) handle treat zero CAVEATS
    % strings as a special case like metadata conventions might require, e.g.
    % replace an empty column CA representing CAVEATS with a scalar CA
    % containing the string "none". This functionality should be handled by
    % bicas.ga.get_output_dataset_GAs() instead.
    %
    function gaCaveats = QRCB_arrays_to_GA_CAVEATS(Qrcbm, Qrcsm)

      assert(isa(Qrcbm, "bicas.proc.QrcbMap"))
      assert(isa(Qrcsm, "bicas.proc.QrcSettingsMap"))
      assert(isequal(Qrcbm.qrcidAr, Qrcsm.qrcidAr))

      gaCaveats = string.empty(0, 1);

      for qrcid = Qrcbm.qrcidAr'
        qrcbAr = Qrcbm.get(qrcid);
        Qrcs   = Qrcsm.get(qrcid);

        if any(qrcbAr)
          % NOTE: Concatenates two column arrays (not column array + scalar).
          gaCaveats = [gaCaveats; Qrcs.gaCaveats];
        end
      end

      % NOTE: CAVEATS should preferably be sorted to yield predictable array of
      % CAVEATS strings.
      gaCaveats = sort(gaCaveats);
    end



    % Convert a (scalar) LxQBM to an array of indices of where the bits are set.
    %
    % RETURN VALUE
    % ============
    % bitPosAr
    %       Column array. Variable-length. Element value 0=LSB.
    %
    function bitPosAr = LxQBM_to_bit_positions(lxqbm)
      % Not vectorized!
      % TODO-DEC: How handle fill value?!
      %   PROPOSAL: Does not need to handle since only intended for analyzing
      %             QRCSs.

      assert(isscalar(lxqbm) & isa(lxqbm, "uint16"))

      bitPosAr = zeros(0, 1);
      for i = 0:15
        if bitget(lxqbm, i+1)
          bitPosAr(end+1, 1) = i;
        end
      end
    end



    % Sets QRCBs to false for saturation QRCs which should not be used,
    % depending on saturation quality scheme.
    %
    % NOTE: QRCBs for both saturation quality schemes are always present
    % regardless of which one is selected.
    %
    % ARGUMENTS
    % =========
    % saturationQualitySchemeId
    %       The saturation quality scheme for which QRCBs should be kept.
    %
    % RETURN VALUE
    % ============
    % Modified copy of argument.
    %
    function Qrcbm = filter_saturation_QRCBs(Qrcbm, saturationQualitySchemeId)
      assert(isstring(saturationQualitySchemeId))
      assert(isa(Qrcbm, "bicas.proc.QrcbMap"))
      assert(isa(Qrcbm, 'handle'))

      switch(saturationQualitySchemeId)
        case "GLOBAL_SATURATION"
          qrcidFalseAr = bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR;

        case "CHANNEL_SATURATION"
          qrcidFalseAr = bicas.const.qrc.Q.GLOBAL_SATURATION_QRCID_AR;

        otherwise
          error("BICAS:ConfigurationBug", ...
            "Illegal argument saturationQualitySchemeId=""%s"".", ...
            saturationQualitySchemeId)
      end

      % Create modified copy of Qrcbm.
      Qrcbm2      = bicas.proc.QrcbMap(Qrcbm.nRecords);
      qrcbFalseAr = false(Qrcbm.nRecords, 1);
      for qrcid = Qrcbm.qrcidAr'
        if ismember(qrcid, qrcidFalseAr)
          qrcbAr = qrcbFalseAr;
        else
          qrcbAr = Qrcbm.get(qrcid);
        end

        Qrcbm2.add(qrcid, qrcbAr)
      end

      % Rename variable.
      Qrcbm = Qrcbm2;
    end



  end    % methods(Static)



end
