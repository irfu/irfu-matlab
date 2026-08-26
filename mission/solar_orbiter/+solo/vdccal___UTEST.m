%
% matlab.unittest automatic test code for solo.vdccal(), i.e. without any kind
% of mocking of calibration data. The tests can therefore not be too
% sophisticated without assuming things about the exact calibration data.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef vdccal___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Unclear how useful this function is and could be. In practice, calibration
    % data and solo.vdccal() input data should to some respect be consistent,
    % since calibration data is effectively based on the input data.
    %
    function test_random_times(T)
      DT_GLOBAL_BEGIN = datetime("2020-02-10T01:02:03Z", "TimeZone", "UTCLeapSeconds");
      DT_GLOBAL_END   = datetime("2026-01-01T01:02:03Z", "TimeZone", "UTCLeapSeconds");
      % NOTE: Execution time depends a lot on this value.
      N_DAYS_JUMP     = 150;

      % ===================================
      % Derive value for some "random" days
      % ===================================
      for DtBegin = DT_GLOBAL_BEGIN : days(N_DAYS_JUMP) : DT_GLOBAL_END
        % Add/subtract some time to not always use UTC midnight.
        DtEnd = DtBegin + hours(1);

        T.test_data(DtBegin, DtEnd)
      end
    end



    % Test dates outside of the calibration data time range.
    function test_out_of_time_range(T)
      T.test_data_begin_dur('2000-01-01T01:02:03Z', hours(1))
      T.test_data_begin_dur('2100-01-01T01:02:03Z', hours(1))
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Call solo.vdccal() with specified input data.
    function [act_DCE_SRF_out, act_PSP_out, act_ScPot_out] = ...
        test(T, VdcTs, EdcTs)
      % ======================
      % CALL CODE TO BE TESTED
      % ======================
      [act_DCE_SRF_out, act_PSP_out, act_ScPot_out, actCodeVerStr, actMatVerStr] = ...
        solo.vdccal(VdcTs, EdcTs, []);

      % NOTE: Using assertion function from BICAS proper to avoid dulicating
      % code.
      bicas.proc.L2L3.ext.assert_vdccal_return_values( ...
        tt2000    =act_DCE_SRF_out.time.ttns, ...
        EdcSrfTs  =act_DCE_SRF_out, ...
        PspTs     =act_PSP_out, ...
        ScpotTs   =act_ScPot_out, ...
        matVerStr =actMatVerStr, ...
        codeVerStr=actCodeVerStr ...
        )
    end



    % Call solo.vdccal() with created data for specified time interval.
    function test_data(T, DtBegin, DtEnd)
      % NOTE: Execution time depends a lot on the default sampling rate.
      DEFAULT_SAMPLING_RATE_HZ = 2;

      tt2000Begin = convertTo(DtBegin, "tt2000");
      tt2000End   = convertTo(DtEnd,   "tt2000");
      tt2000Ar    = [tt2000Begin : int64(1e9/DEFAULT_SAMPLING_RATE_HZ) : tt2000End]';

      nTimestamps = numel(tt2000Ar);

      % ====================================
      % TEST 1: Nonsense data, always finite
      % ====================================
      dataAr      = zeros(nTimestamps, 3);
      VdcTs       = irf.ts_vec_xyz(EpochTT(tt2000Ar), dataAr);
      EdcTs       = irf.ts_vec_xyz(EpochTT(tt2000Ar), dataAr);
      T.test(VdcTs, EdcTs);

      % ==========================
      % TEST 2: Nonsense data, NaN
      % ==========================
      dataAr      = zeros(nTimestamps, 3) * NaN;
      VdcTs       = irf.ts_vec_xyz(EpochTT(tt2000Ar), dataAr);
      EdcTs       = irf.ts_vec_xyz(EpochTT(tt2000Ar), dataAr);
      [act_DCE_SRF_out, act_PSP_out, act_ScPot_out] = T.test(VdcTs, EdcTs);
      assert(all(isnan(act_DCE_SRF_out.data), 'all'))
      assert(all(isnan(act_PSP_out.data    ), 'all'))
      assert(all(isnan(act_ScPot_out.data  ), 'all'))

      % ======================================================
      % TEST 3: Nonsense data, partially NaN, partially finite
      % ======================================================
      % NOTE: Assumes sufficiently many timestamps for the fake timestamps.
      dataAr = [...
        NaN,   2,   3; ...
        1,   NaN,   3; ...
        1,     2, NaN; ...
        NaN, NaN,   3; ...
        1,   NaN, NaN; ...
        NaN,   2, NaN; ...
        ];
      VdcTs = irf.ts_vec_xyz(EpochTT(tt2000Ar(1:size(dataAr, 1))), dataAr);
      EdcTs = irf.ts_vec_xyz(EpochTT(tt2000Ar(1:size(dataAr, 1))), dataAr);
      T.test(VdcTs, EdcTs);
    end



    % Call solo.vdccal() with created data for specified time interval.
    % Primarily intended for hardcoded timestamps.
    function test_data_begin_dur(T, dateStrBegin, Duration)
      assert(isa(Duration, "duration"))

      DtBegin = datetime(dateStrBegin, "TimeZone", "UTCLeapSeconds");
      DtEnd   = DtBegin + Duration;

      T.test_data(DtBegin, DtEnd)
    end



  end    % methods(Access=private)



end
