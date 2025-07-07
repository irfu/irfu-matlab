%
% Collection of reusable, generic code relating to setting quality ZVs.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual
  % PROPOSAL: Class for containers.Map QRCID-->QRCB array.
  %   PRO: Simplifies documentation (description of arguments).
  %   PRO: Easier assertions.
  %     PRO: Shorter assertions when used in code. Reused assertions in class.
  %   PRO: Natural grouping of related functionality including tests.
  %   PRO: Simplifies tests
  %     PRO: Easier to initialize empty map.
  %     CON: Must have tests for equality.
  %   CON: Must implement support for equality for test code.
  %     CON: Easy to implement equality. Can just delegate to isequal(map1,
  %          map2).
  %   NOTE: ~Needs abbreviation?
  %     PROPOSAL: QRCBM
  %   NOTE: ~Must be handle class.
  %   TODO-NI: Class name?
  %     ~QRCID, ~QRCB, ~map
  %     ~arrays
  %     QrcbMap
  %       QRCBM=
  %     QrcBitMap
  %       CON: Sounds too much like bitmap (image).
  %     QrcbArraysMap
  %       QAM=



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % IMPLEMENTATION NOTE: Without allQrcidCa, the function can not create a
    % return value map that contains keys for all QRCIDs, in case the
    % NsoTable does not contain all QRCIDs.
    %
    % ARGUMENTS
    % =========
    % requestedQrcidAr
    %       Column array of QRCIDs for which the return value shall contain
    %       data. This is a way of filtering the content of the NSO table. The
    %       function will return default values (QRCB=false) for QRCIDs for
    %       which there are no events which overlap with the specified time
    %       stamps.
    %
    %
    % RETURN VALUE
    % ============
    % QrcbMap
    %       containers.Map. QRCID->logical array
    %       Contains keys for all QRCIDs specified in allQrcidCa, not just
    %       those present in NsoTable.
    %
    function QrcbMap = NSO_table_to_QRCB_arrays(...
        requestedQrcidAr, NsoTable, tt2000Ar, L)

      % PROPOSAL: Redefine allQrcidAr to only be those QRCIDs for which one
      %           wants to be used for populating the return value QrcbMap.
      %   PRO: Useful for separating relevant QRCs, e.g. for only L2, only L3,
      %        only L3 density etc.
      %   CON: Can not use allQrcidAr for asserting only legal QRCIDs in NSO
      %        table.
      %   CON-PROPOSAL: Separate function for filtering QrcbMap wrt. QRCIDs.

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
      % Initialize "empty" QrcbMap (elements=false) for all requested QRCIDs
      %----------------------------------------------------------------------
      % IMPLEMENTATION NOTE: valueType=logical implies scalar (sic!) and can
      %                      therefore not be used.
      QrcbMap = containers.Map('keyType', 'char', 'valueType', 'any');
      for i = 1:numel(requestedQrcidAr)
        QrcbMap(requestedQrcidAr(i)) = false(size(tt2000Ar));
      end

      %-----------------------------------------------------------------------
      % Iterate over NSO events and set QRCBs for resp. QRCIDs and timestamps
      %-----------------------------------------------------------------------
      for kLe = 1:nLe
        % Index into GLOBAL NSO events table.
        iGe        = NsoEventMatchAr(kLe).iNsoEvent;
        eventQrcid = NsoEventMatchAr(kLe).qrcid;
        % Logical indices into tt2000Ar.
        qrbcAr     = NsoEventMatchAr(kLe).qrcbAr;

        %=====================================================================
        % Log relevant NSO events by referring to the GLOBAL NSO events table
        %=====================================================================
        L.logf('info', '    %s -- %s %s', ...
          bicas.utils.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGe), 9), ...
          bicas.utils.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array( iGe), 9), ...
          eventQrcid);

        % ASSERTION
        % NOTE: Not perfect assertion on legal QRCIDs since the code only checks
        % those relevant for the data (time interval) currently processed.
        % (Therefore also checks all QRCIDs when reads NSO table.)
        assert(ismember(eventQrcid, requestedQrcidAr), ...
          'Can not interpret QRCID "%s".', eventQrcid)

        %======================================
        % Set corresponding QRC array elements
        %======================================
        QrcbMap(eventQrcid) = QrcbMap(eventQrcid) | qrbcAr;
      end    % for
    end



    % Given an array(s) of QRC bits (one bit per CDF record and QRC), translate
    % that into
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
    % nRec
    %       Number of CDF records (rows).
    %       IMPLEMENTATION NOTE: Needed for handling the case of zero QRCIDs.
    % QrcbMap
    %       containers.Map: QRCID-->Logical column array
    %       Array element is set when the corresponding QRC applies.
    % QrcsMap
    %       containers.Map: QRCID-->bicas.proc.QrcSetting
    %
    %
    % RETURN VALUES
    % =============
    % QUALITY_FLAG
    %       QUALITY_FLAG max value wrt. to QRCs handled in this function.
    % Lx_QUALITY_BITMASK
    %       L*_QUALITY_BITMASK value wrt. to QRCs handled in this function.
    %       Refers to L2_QUALITY_BITMASK or L3_QUALITY_BITMASK depending on
    %       context.
    %
    function [QUALITY_FLAG, Lx_QUALITY_BITMASK] = QRCB_arrays_to_quality_ZVs(...
        nRec, QrcbMap, QrcsMap)

      % PROPOSAL: Keys in QrcbMap and QrcsMap do not have to be identical.
      %   -- IMPLEMENTED
      %   CON: Loses fail-check.
      %   PRO: Easier to call.
      %     PRO: Can just submit ~global map of  QRCSs

      %irf.assert.castring_sets_equal(QrcbMap.keys, QrcsMap.keys)

      % Create "empty" quality variable arrays, with max possible quality
      % (QUALITY_FLAG max, quality bits=0), which can then later be "lowered"
      % if necessary.
      QUALITY_FLAG       = ones( nRec, 1, 'uint8' ) * bicas.const.QUALITY_FLAG_MAX;
      Lx_QUALITY_BITMASK = zeros(nRec, 1, 'uint16');

      %=================================
      % Iterate over QRCIDs in argument
      %=================================
      qrcidCa = QrcbMap.keys();
      for i = 1:numel(qrcidCa)
        qrcid  = qrcidCa{i};
        Qrcs   = QrcsMap(qrcid);
        qrcbAr = QrcbMap(qrcid);

        assert(isa(qrcbAr, 'logical'))
        assert(isequal( size(qrcbAr), [nRec, 1] ))

        %------------------
        % Set QUALITY_FLAG
        %------------------
        % IMPLEMENTATION NOTE: Only adjusts relevant indices since the
        % operation is more natural (simpler) that way.
        QUALITY_FLAG(qrcbAr) = min(...
          QUALITY_FLAG(qrcbAr), ...
          Qrcs.QUALITY_FLAG);

        %------------------------
        % Set Lx_QUALITY_BITMASK
        %------------------------
        Lx_QUALITY_BITMASK = bitor(...
          Lx_QUALITY_BITMASK, ...
          Qrcs.Lx_QUALITY_BITMASK * uint16(qrcbAr));
      end
    end



  end    % methods(Static)



end
