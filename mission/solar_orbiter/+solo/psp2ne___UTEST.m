%
% matlab.unittest automatic test code for solo.psp2ne(), i.e. without any kind
% of mocking of calibration data. The tests can therefore not be too
% sophisticated without assuming things about the exact calibration data.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef psp2ne___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  % Additional properties of testCase objects. Needed for setup and teardown
  % methods which store/read their own data from the testCase object.
  properties
    L
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(T)
      T.L = bicas.Logger('HUMAN_READABLE', false);
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % NOTE: Some test data relies on the hardcoded calibration data in
    % solo.psp2ne() and may have to be updated when solo.psp2ne() is updated.
    %
    % NOTE: Does not currently test the small gaps between time intervals in
    % hardcoded psp2ne() calibration data.
    %
    % NOTE: Some tests may in principle be better with higher sampling rate
    % (checking 1-second gaps between hardcoded calibration data time
    % intervals).
    %
    % NOTE: Return value validation should preferably be consistent with tests
    % in bicas.proc.L2L3.ext where BICAS nominally calls solo.psp2ne().
    %
    function test_0(T)
      DT_GLOBAL_BEGIN = datetime("2020-02-10T00:00:00Z", "TimeZone", "UTCLeapSeconds");
      DT_GLOBAL_END   = datetime("2026-01-01T00:00:00Z", "TimeZone", "UTCLeapSeconds");
      % NOTE: Execution time depends a lot on this value.
      N_DAYS_JUMP     = 150;

      % ===================================
      % Derive value for some "random" days
      % ===================================
      %tic
      for DtBegin = DT_GLOBAL_BEGIN : days(N_DAYS_JUMP) : DT_GLOBAL_END
        % Add/subtract some time to not always use UTC midnight.
        WIGGLE_TIME = minutes(1)*2;
        DtEnd = DtBegin + hours(1);

        T.test(DtBegin - WIGGLE_TIME, DtEnd, "MIXED")
        T.test(DtBegin + WIGGLE_TIME, DtEnd, "MIXED")
      end
      %toc



      % ===============================
      % Test specific dates/time ranges
      % ===============================

      % Call for entire mission.
      T.test(DT_GLOBAL_BEGIN, DT_GLOBAL_END, "MIXED", 2/86400)

      % Test dates which are definitively outside of hardcoded calibration data
      T.test_begin_dur("2020-01-01T00:00:00Z", hours(1), "NAN")
      T.test_begin_dur("2120-01-01T00:00:00Z", hours(1), "NAN")

      % Test AddEntry() with real numbers
      % ---------------------------------
      % AddEntry('2021-01-01T00:00:00Z/2021-01-06T23:59:59Z',[0.5905  4.0923]);
      % AddEntry('2021-01-07T00:00:00Z/2021-01-12T05:49:59Z',[0.6730  4.1837]); %2
      % AddEntry('2021-01-12T05:50:00Z/2021-01-17T23:59:59Z',[0.7462  4.5630]); %3
      T.test_begin_dur("2021-01-06T00:00:00Z", days(1),  "FINITE")
      T.test_begin_dur("2021-01-12T00:00:00Z", days(1),  "FINITE")

      % Test AddEntry() with complex numbers/2-fold approximations
      % ----------------------------------------------------------
      % AddEntry('2021-02-17T00:00:00Z/2021-03-20T01:29:59Z',[0.7651  3.4971]); %5
      % AddEntry('2021-03-20T01:30:00Z/2021-03-22T19:29:59Z',... %6
      %   [0.6460 + 3.5047i   0.3899 + 3.7683i],1.0297);
      % AddEntry('2021-03-22T19:30:00Z/2021-04-04T03:59:59Z',[0.7884  3.3714]); %7
      T.test_begin_dur("2021-03-20T00:00:00Z", hours(1), "FINITE")
      T.test_begin_dur("2021-03-21T00:00:00Z", hours(1), "FINITE")
      T.test_begin_dur("2021-03-22T00:00:00Z", hours(1), "FINITE")

      % Test interior period explicitly hardcoded to be NaN
      % ---------------------------------------------------
      % AddEntry('2023-09-08T00:00:00Z/2025-02-28T23:59:59Z',[NaN,    NaN   ]);    % NOTE: ~18 MONTHS calibration data gap!!
      T.test_begin_dur("2024-01-01T00:30:00Z", hours(1), "NAN")
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Call solo.psp2ne() with created data for specified time interval.
    % Primarily intended for generated timestamps.
    %
    % ARGUMENTS
    % =========
    % ExpRvId
    %       String constant. Specifies what kind of expected return values.
    %
    function test(T, DtBegin, DtEnd, ExpRvId, varargin)

      % NOTE: Execution time depends a lot on the default sampling rate.
      DEFAULT_SAMPLING_RATE_HZ = 2;

      % NOTE: Optional "keyword argument" "samplingRateHz" does not seem to be
      %       used.
      p = inputParser;
      p.addOptional("samplingRateHz", DEFAULT_SAMPLING_RATE_HZ)
      p.parse(varargin{:})
      samplingRateHz = p.Results.samplingRateHz;

      assert(DtBegin < DtEnd)
      assert(isstring(ExpRvId))



      tt2000Begin = convertTo(DtBegin, "tt2000");
      tt2000End   = convertTo(DtEnd,   "tt2000");

      tt2000Ar    = [tt2000Begin : int64(1e9/samplingRateHz) : tt2000End]';
      nTimestamps = numel(tt2000Ar);
      dataAr      = zeros(nTimestamps, 1);

      PspTs       = TSeries(EpochTT(tt2000Ar), dataAr);



      % CALL CODE TO BE TESTED
      [ActNeScpTs, ActNeScpQualityBitTs, actCodeVerStr] = solo.psp2ne(PspTs);



      % NOTE: Using assertion function from BICAS proper to avoid dulicating
      % code.
      bicas.proc.L2L3.ext.assert_psp2ne_return_values( ...
        PspTs.time.ttns, ActNeScpTs, ActNeScpQualityBitTs, actCodeVerStr, T.L)

      switch(ExpRvId)
        case "NAN"
          assert(all(isnan(ActNeScpTs.data)))

          % NOTE: Condition reflects actual behaviour. Not sure if this is
          % desirable behaviour though.
          assert(all(ActNeScpQualityBitTs.data == 0))

        case "FINITE"
          assert(all(isfinite(ActNeScpTs.data)))

        case "MIXED"
          ;   % Do nothing

        otherwise
          error("Illegal ExpId")
      end
    end



    % Call solo.psp2ne() with created data for specified time interval.
    % Primarily intended for hardcoded timestamps.
    function test_begin_dur(T, dateStrBegin, Duration, varargin)
      assert(isa(Duration, "duration"))

      DtBegin = datetime(dateStrBegin, "TimeZone", "UTCLeapSeconds");
      DtEnd   = DtBegin + Duration;

      T.test(DtBegin, DtEnd, varargin{:})
    end



  end    % methods(Access=private)



end
