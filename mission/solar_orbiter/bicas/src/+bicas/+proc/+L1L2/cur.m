%
% Code related to bias currents, broken out from bicas.proc.L1L2.dc.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef cur



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Calibrate currents on the format they are found in CDFs.
    %
    function currentAampere = calibrate_bias_currents( ...
        curTt2000Ar, currentSampere, sciTt2000Ar, Ccal, Bso, L)

      assert(isa(Ccal, "bicas.proc.L1L2.cal.CurrentCalibrationAbstract"))

      currentSampere = bicas.proc.L1L2.cur.convert_CUR_to_CUR_on_SCI_TIME(...
        curTt2000Ar, currentSampere, sciTt2000Ar, Bso, L);
      currentTm      = Ccal.calibrate_current_sampere_to_TM(currentSampere);

      currentAampere = nan(size(currentSampere));    % Preallocate.
      iCalibL        = Ccal.get_BIAS_calibration_time_index_L(sciTt2000Ar);

      L.logf('info', ...
        ['Calibrating currents -', ...
        ' One sequence of records with identical settings at a time.'])
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(iCalibL);
      for iSs = 1:nSs
        iRec1    = iRec1Ar(iSs);
        iRec2    = iRec2Ar(iSs);
        iRecords = iRec1:iRec2;

        L.logf('info', 'Records %8i-%8i : %s -- %s', ...
          iRec1, iRec2, ...
          bicas.utils.TT2000_to_UTC_str(sciTt2000Ar(iRec1), 9), ...
          bicas.utils.TT2000_to_UTC_str(sciTt2000Ar(iRec2), 9))

        for iAnt = 1:3
          %--------------------
          % CALIBRATE CURRENTS
          %--------------------
          currentAampere(iRecords, iAnt) = Ccal.calibrate_current_TM_to_aampere(...
            currentTm(iRecords, iAnt), iAnt, iCalibL(iRecords));
        end
      end    % for

    end    % function



    function currentSampere = convert_CUR_to_CUR_on_SCI_TIME(...
        curTt2000Ar, currentSampereAr, sciTt2000Ar, Bso, L)

      % PROPOSAL: Better name.
      %   PRO: Does not imply unit, or changing units (set nanoampere-->set
      %        ampere).
      %     CON: Argument and return value name implies unit.

      % ASSERTIONS
      irf.assert.sizes(...
        curTt2000Ar,      [-1], ...
        currentSampereAr, [-1, 3])
      bicas.utils.assert_ZV_Epoch(curTt2000Ar)
      bicas.utils.assert_ZV_Epoch(sciTt2000Ar)



      %===================================================================
      % CDF ASSERTION: CURRENT data begins before SCI data (i.e. there is
      % enough CURRENT data).
      %===================================================================
      if ~(min(curTt2000Ar) <= min(sciTt2000Ar))
        curRelativeSec    = 1e-9 * (min(curTt2000Ar) - min(sciTt2000Ar));
        sciEpochUtcStr    = bicas.utils.TT2000_to_UTC_str(min(sciTt2000Ar), 9);
        curEpochMinUtcStr = bicas.utils.TT2000_to_UTC_str(min(curTt2000Ar), 9);

        [settingValue, settingKey] = Bso.get_fv(...
          'PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY');

        anomalyDescrMsg = sprintf(...
          ['Bias current data begins %g s (%s) AFTER voltage data begins (%s).', ....
          ' Can therefore not determine currents for all voltage timestamps.'], ...
          curRelativeSec, curEpochMinUtcStr, sciEpochUtcStr);

        bicas.default_anomaly_handling(L, settingValue, settingKey, ...
          'ERROR_WARNING_ILLEGAL_SETTING', ...
          anomalyDescrMsg, 'BICAS:SWMProcessing')
      end



      %====================================================================
      % CDF ASSERTION: Epoch increases (not monotonically)
      % --------------------------------------------------
      % NOTE: bicas.proc.L1L2.cur.zv_TC_to_current() checks (and handles)
      % that Epoch increases monotonically, but only for each antenna
      % separately (which does not capture all cases). Therefore checks
      % that Epoch is (non-monotonically) increasing.
      % Ex: Timestamps, iAntenna = mod(iRecord,3): 1,2,3,5,4,6
      %       ==> Monotonically increasing sequences for each antenna
      %           separately, but not even increasing when combined.
      %====================================================================
      assert(issorted(curTt2000Ar), ...
        'BICAS:DatasetFormat', ...
        'CURRENT timestamps zVar Epoch does not increase (all antennas combined).')

      % NOTE: bicas.proc.L1L2.cur.zv_TC_to_current() checks that Epoch
      % increases monotonically.
      currentNanoSampere = [];
      for iAnt = 1:3
        currentNanoSampere(:,iAnt) = bicas.proc.L1L2.cur.zv_TC_to_current(...
          curTt2000Ar, currentSampereAr(:,iAnt), sciTt2000Ar, L, Bso);
      end

      currentSampere = 1e-9 * currentNanoSampere;
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Wrapper around solo.hwzv.CURRENT_ZV_to_current_interpolate() for anomaly
    % handling.
    function sciZv_IBIASx = zv_TC_to_current(...
        curZv_Epoch, curZv_IBIAS_x, sciZv_Epoch, L, Bso)

      %====================
      % Calibrate currents
      %====================
      [sciZv_IBIASx, duplicateAnomaly] = ...
        solo.hwzv.CURRENT_ZV_to_current_interpolate(...
        curZv_Epoch, ...
        curZv_IBIAS_x, ...
        sciZv_Epoch);



      if duplicateAnomaly
        %====================================================
        % Handle anomaly: Non-monotonically increasing Epoch
        %====================================================
        [settingValue, settingKey] = Bso.get_fv(...
          'INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY');
        anomalyDescriptionMsg = [...
          'Bias current data contain duplicate settings, with', ...
          ' identical timestamps', ...
          ' and identical bias settings on the same antenna.'];

        switch(settingValue)
          case 'REMOVE_DUPLICATES'
            bicas.default_anomaly_handling(L, ...
              settingValue, settingKey, 'OTHER', ...
              anomalyDescriptionMsg)
            L.log('warning', ...
              ['Removed duplicated bias current settings with', ...
              ' identical timestamps on the same antenna.'])

          otherwise
            bicas.default_anomaly_handling(L, ...
              settingValue, settingKey, 'ERROR_ILLEGAL_SETTING', ...
              anomalyDescriptionMsg, ...
              'BICAS:SWMProcessing:DatasetFormat')
        end
      end

    end    % bicas.proc.L1L2.cur.zv_TC_to_current



  end    % methods(Static, Access=private)



end
