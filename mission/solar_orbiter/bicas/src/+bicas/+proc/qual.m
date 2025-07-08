%
% Collection of reusable, generic code relating to setting quality ZVs.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % IMPLEMENTATION NOTE: Without requestedQrcidAr, the function can not create
    % a return value map that contains keys for all QRCIDs, in case the NsoTable
    % does not contain all QRCIDs.
    %
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
    %       Contains keys for all QRCIDs specified in requestedQrcidAr, not just
    %       those present in NsoTable.
    %
    function QrcbMap = NSO_table_to_QRCB_map(...
        requestedQrcidAr, NsoTable, tt2000Ar, L)

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
      QrcbMap = bicas.proc.QrcbMap(numel(tt2000Ar));
      for i = 1:numel(requestedQrcidAr)
        QrcbMap.add(requestedQrcidAr(i), false(size(tt2000Ar)));
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
          QrcbMap.set(eventQrcid, QrcbMap.get(eventQrcid) | qrbcAr);
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
    % nRec
    %       Number of CDF records (rows).
    %       IMPLEMENTATION NOTE: Needed for handling the case of zero QRCIDs.
    % QrcbMap
    %       containers.Map: QRCID-->Logical column array
    %       Array element is set when the corresponding QRC applies.
    % QrcsMap
    %       containers.Map: QRCID-->bicas.proc.QrcSetting
    %       Must contain at least the QRCIDs in QrcbMap.
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

      assert(isa(QrcbMap, "bicas.proc.QrcbMap"))

      % Create "empty" quality variable arrays, with max possible quality
      % (QUALITY_FLAG max, Lx_QUALITY_BITMASK=0), which can then later be
      % "lowered" if necessary.
      QUALITY_FLAG       = ones( nRec, 1, 'uint8' ) * bicas.const.QUALITY_FLAG_MAX;
      Lx_QUALITY_BITMASK = zeros(nRec, 1, 'uint16');

      %=================================
      % Iterate over QRCIDs in argument
      %=================================
      for qrcid = QrcbMap.qrcidAr'
        Qrcs   = QrcsMap(qrcid);
        qrcbAr = QrcbMap.get(qrcid);

        assert(isa(qrcbAr, 'logical'))
        assert(isequal( size(qrcbAr), [nRec, 1] ))

        %------------------
        % Set QUALITY_FLAG
        %------------------
        % IMPLEMENTATION NOTE: Only adjusts relevant indices since the operation
        % is more natural (simpler, shorter) that way. min() should only be
        % applied to the indices where QRCB=true.
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
