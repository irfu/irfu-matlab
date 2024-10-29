%
% Code related to bias currents, broken out from bicas.proc.L1L2.dc.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef cur
  % PROPOSAL: Automatic test code.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function currentAAmpere = calibrate_bias_currents(sciEpoch, InCurPd, Cal, Bso, L)
      currentSAmpere = bicas.proc.L1L2.cur.convert_CUR_to_CUR_on_SCI_TIME(...
        sciEpoch, InCurPd, Bso, L);
      currentTm      = bicas.proc.L1L2.cal.Cal.calibrate_current_sampere_to_TM(currentSAmpere);

      currentAAmpere          = nan(size(currentSAmpere));    % Preallocate.
      iCalibLZv               = Cal.get_BIAS_calibration_time_L(sciEpoch);
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(iCalibLZv);
      L.logf('info', ...
        ['Calibrating currents -', ...
        ' One sequence of records with identical settings at a time.'])
      for iSs = 1:nSs
        iRec1    = iRec1Ar(iSs);
        iRec2    = iRec2Ar(iSs);
        iRecords = iRec1:iRec2;

        L.logf('info', 'Records %8i-%8i : %s -- %s', ...
          iRec1, iRec2, ...
          bicas.utils.TT2000_to_UTC_str(sciEpoch(iRec1), 9), ...
          bicas.utils.TT2000_to_UTC_str(sciEpoch(iRec2), 9))

        for iAnt = 1:3
          %--------------------
          % CALIBRATE CURRENTS
          %--------------------
          currentAAmpere(iRecords, iAnt) = Cal.calibrate_current_TM_to_aampere(...
            currentTm(iRecords, iAnt), iAnt, iCalibLZv(iRecords));
        end
      end    % for

    end    % function



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function currentSAmpere = convert_CUR_to_CUR_on_SCI_TIME(...
        sciEpoch, InCur, Bso, L)

      % PROPOSAL: Change function name. process_* implies converting struct-->struct.

      % ASSERTIONS
      assert(isa(InCur, 'bicas.InputDataset'))



      %===================================================================
      % CDF ASSERTION: CURRENT data begins before SCI data (i.e. there is
      % enough CURRENT data).
      %===================================================================
      if ~(min(InCur.Zv.Epoch) <= min(sciEpoch))
        curRelativeSec    = 1e-9 * (min(InCur.Zv.Epoch) - min(sciEpoch));
        sciEpochUtcStr    = bicas.utils.TT2000_to_UTC_str(min(sciEpoch),       9);
        curEpochMinUtcStr = bicas.utils.TT2000_to_UTC_str(min(InCur.Zv.Epoch), 9);

        [settingValue, settingKey] = Bso.get_fv(...
          'PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY');

        anomalyDescrMsg = sprintf(...
          ['Bias current data begins %g s (%s) AFTER voltage data begins (%s).', ....
          ' Can therefore not determine currents for all voltage timestamps.'], ...
          curRelativeSec, curEpochMinUtcStr, sciEpochUtcStr);

        bicas.default_anomaly_handling(L, settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
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
      assert(issorted(InCur.Zv.Epoch), ...
        'BICAS:DatasetFormat', ...
        'CURRENT timestamps zVar Epoch does not increase (all antennas combined).')

      % NOTE: bicas.proc.L1L2.cur.zv_TC_to_current() checks that Epoch
      % increases monotonically.
      currentNanoSAmpere = [];
      currentNanoSAmpere(:,1) = bicas.proc.L1L2.cur.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_1, sciEpoch, L, Bso);
      currentNanoSAmpere(:,2) = bicas.proc.L1L2.cur.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_2, sciEpoch, L, Bso);
      currentNanoSAmpere(:,3) = bicas.proc.L1L2.cur.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_3, sciEpoch, L, Bso);

      currentSAmpere = 1e-9 * currentNanoSAmpere;
    end



    % Wrapper around solo.hwzv.CURRENT_ZV_to_current_interpolate for
    % anomaly handling.
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
