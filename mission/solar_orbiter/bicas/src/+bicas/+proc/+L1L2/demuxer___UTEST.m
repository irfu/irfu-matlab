%
% matlab.unittest automatic test code for bicas.proc.L1L2.demuxer.
%
% Could be improved but unsure how much is meaningful. Seems to complicated.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08, using older test code.
%
classdef demuxer___UTEST < matlab.unittest.TestCase



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function AsidTestSamplesSrm = create_channel_test_data(sampleSize)
      A = bicas.proc.L1L2.const.C.ASID_DICT;

      % (iChannel, 1) = ASID.
      % (iChannel, 2) = Sample value (that is consistent with other
      %                 channels).
      TEST_DATA_CA = { ...
        A("DC_V1"),  10; ...
        A("DC_V2"),  11; ...
        A("DC_V3"),  13; ...
        A("DC_V12"), 10-11; ...
        A("DC_V13"), 10-13; ...
        A("DC_V23"), 11-13; ...
        A("AC_V12"), 45-56; ...
        A("AC_V13"), 45-69; ...
        A("AC_V23"), 56-69 ...
        };

      % Multiply the sample values with matrix to test multiple records
      % with "snapshots" (SPR>1).
      TEST_DATA_CA(:, 2) = cellfun( ...
        @(x) (x * ones(sampleSize)), TEST_DATA_CA(:, 2), ...
        'UniformOutput', false);

      AsidTestSamplesSrm = bicas.utils.SameRowsMap( ...
        "uint8", sampleSize(1), 'EMPTY');
      for i = 1:size(TEST_DATA_CA, 1)
        AsidTestSamplesSrm.add( ...
          TEST_DATA_CA{i, 1}, ...
          TEST_DATA_CA{i, 2})
      end
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Test two functions in combination:
    %   bicas.proc.L1L2.demuxer.get_routings()
    %   bicas.proc.L1L2.demuxer.relabel_reconstruct_samples_BLTS_to_ASR_subsequence()
    %
    % IMPLEMENTATION NOTE: The design is for historical reasons before the
    % underlying code (old function) was split up in two functions.
    %
    function test_get_routings_relabel_reconstruct_samples_BLTS_to_ASR_subse(testCase)
      % PROBLEM: Too large test function.
      % PROPOSAL: Separate test code for the two functions.
      %   CON: Combining both ensures that tests only test real combinations
      %        of routings.

      A = bicas.proc.L1L2.const.C.ASID_DICT;
      R = bicas.proc.L1L2.const.C.ROUTING_DICT;

      % =========
      % Test data
      % =========
      TEST_DATA_UNKNOWN = [999];   % Data from unknown source.
      SAMPLES_SIZE      = [3,2];

      nRows = SAMPLES_SIZE(1);
      AsidTestSamplesSrm = testCase.create_channel_test_data(SAMPLES_SIZE);



      % Test any BDM.
      function test(bdmFloatNan, dlrFloatNan, bltsSamplesAVolt, ExpRoutingArray, ExpAsrSamplesAVoltSrm)
        assert(numel(ExpRoutingArray) == 5)

        dlrFpa = bicas.utils.FPArray.floatNan2logical(dlrFloatNan);
        bdmFpa = bicas.utils.FPArray.floatNan2int(bdmFloatNan, 'uint8');

        % CALL FUNCTIONS
        ActRoutingArray       = bicas.proc.L1L2.demuxer.get_routings(...
          bdmFpa, dlrFpa);
        ActAsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.relabel_reconstruct_samples_BLTS_to_ASR_subsequence(...
          [ActRoutingArray.sdid], bltsSamplesAVolt);

        % ASSERTIONS
        testCase.assertEqual(ActRoutingArray, ExpRoutingArray)
        testCase.assertTrue(ActAsrSamplesAVoltSrm == ExpAsrSamplesAVoltSrm)
      end



      % Function for testing BDM 0-4. All those all those map (route) ASR
      % to ASR (no GNS, no 2.5V REF, no "unknown", no "nowhere").
      function test_BDM01234(bdmFloatNan, dlrFloatNan, ExpRoutingArray, ExpAsrSamplesAVoltSrm)
        assert(isa(ExpAsrSamplesAVoltSrm, "bicas.utils.SameRowsMap"))
        assert(ismember(bdmFloatNan, [0:4]))

        % Autogenerate bltsSamplesCa (test argument) using
        % ExpRoutingArray (only possible for BDM 0-4.
        tempBltsSamplesAVolt = gen_BLTS_samples(ExpRoutingArray);

        test(bdmFloatNan, dlrFloatNan, tempBltsSamplesAVolt, ExpRoutingArray, ExpAsrSamplesAVoltSrm)
      end



      function tempBltsSamplesAVolt = gen_BLTS_samples(RoutingArray)
        tempBltsSamplesAVolt = zeros(SAMPLES_SIZE);

        for i = 1:numel(RoutingArray)
          Routing = RoutingArray(i);
          if bicas.proc.L1L2.const.SSID_is_ASR(Routing.ssid)
            tempBltsSamplesAVolt(:, :, i) = AsidTestSamplesSrm(...
              bicas.proc.L1L2.const.SSID_ASR_to_ASID(Routing.ssid));
          else
            tempBltsSamplesAVolt(:, :, i) = TEST_DATA_UNKNOWN;
          end
        end
      end



      % Create samples per ASID using constants. Arguments determine for
      % which ASIDs samples should be NaN instead of data.
      %
      % varargin{i} == 0 or 1. Determines whether constant or NaN will
      % be used.
      function AsSrm = get_ASR_samples(varargin)
        assert(nargin == 9)

        % Define which varargin{i} corresponds to which ASID
        % --------------------------------------------------
        % NOTE: This is not (necessarily) the same as
        % bicas.proc.L1L2.const.C.ASID_DICT.values.
        ARGS_ASID_ARRAY = A([
          "DC_V1",  "DC_V2",  "DC_V3",  ...
          "DC_V12", "DC_V13", "DC_V23", ...
          "AC_V12", "AC_V13", "AC_V23" ...
          ]);
        AsSrm = bicas.utils.SameRowsMap("uint8", nRows, 'EMPTY');

        for iAsid = 1:9
          asid = ARGS_ASID_ARRAY(iAsid);

          samplesAVolt = AsidTestSamplesSrm(asid);
          if ~varargin{iAsid}
            samplesAVolt = samplesAVolt * NaN;
          end
          AsSrm.add(asid, samplesAVolt);
        end
      end



      % ==========================
      % bdm = 0, dlr = [0, 1, NaN]
      % ==========================
      test_BDM01234(...
        0, 0, ...
        R(["DC_V1", "DC_V12", "DC_V23", "AC_V12", "AC_V23"]), ...
        get_ASR_samples(1,1,1, 1,1,1, 1,1,1)...
        )
      test_BDM01234(...
        0, 1, ...
        R(["DC_V1", "DC_V13", "DC_V23", "AC_V13", "AC_V23"]), ...
        get_ASR_samples(1,1,1, 1,1,1, 1,1,1)...
        )
      test_BDM01234(...
        0, NaN, ...
        R(["DC_V1", "UNKNOWN_TO_NOWHERE", "DC_V23", "UNKNOWN_TO_NOWHERE", "AC_V23"]), ...
        get_ASR_samples(1,0,0, 0,0,1, 0,0,1)...
        )

      % =======================
      % bdm = 1, dlr = [0, NaN]
      % =======================
      test_BDM01234(...
        1, 1, ...
        R(["DC_V2", "DC_V3", "DC_V23", "AC_V13", "AC_V23"]), ...
        get_ASR_samples(0,1,1, 0,0,1, 1,1,1) ...
        )
      test_BDM01234(...
        1, NaN, ...
        R(["DC_V2", "DC_V3", "DC_V23", "UNKNOWN_TO_NOWHERE", "AC_V23"]), ...
        get_ASR_samples(0,1,1, 0,0,1, 0,0,1) ...
        )

      % ==============================
      % bdm = 4 (calibration), dlr = 1
      % ==============================
      test_BDM01234(...
        4, 0, ...
        R(["DC_V1", "DC_V2", "DC_V3", "AC_V12", "AC_V23"]), ...
        get_ASR_samples(1,1,1, 1,1,1, 1,1,1) ...
        )



      % ==============
      % BDM = 5, DLR 1
      % ==============
      bltsSamplesAVolt(:, :, 1) = AsidTestSamplesSrm(A("DC_V1"));
      bltsSamplesAVolt(:, :, 2) = AsidTestSamplesSrm(A("DC_V2"));
      bltsSamplesAVolt(:, :, 3) = AsidTestSamplesSrm(A("DC_V3"));
      bltsSamplesAVolt(:, :, 4) = AsidTestSamplesSrm(A("AC_V13"));
      bltsSamplesAVolt(:, :, 5) = AsidTestSamplesSrm(A("AC_V23"));
      test(5, 1, ...
        bltsSamplesAVolt, ...
        R(["REF25V_TO_DC_V1", "REF25V_TO_DC_V2", "REF25V_TO_DC_V3", "AC_V13", "AC_V23"]), ...
        get_ASR_samples(1,1,1, 1,1,1, 1,1,1))

      % ==============
      % BDM = 6, DLR 0
      % ==============
      bltsSamplesAVolt(:, :, 1) = AsidTestSamplesSrm(A("DC_V1"));
      bltsSamplesAVolt(:, :, 2) = AsidTestSamplesSrm(A("DC_V2"));
      bltsSamplesAVolt(:, :, 3) = AsidTestSamplesSrm(A("DC_V3"));
      bltsSamplesAVolt(:, :, 4) = AsidTestSamplesSrm(A("AC_V12"));
      bltsSamplesAVolt(:, :, 5) = AsidTestSamplesSrm(A("AC_V23"));
      test(6, 0, ...
        bltsSamplesAVolt, ...
        R(["GND_TO_DC_V1", "GND_TO_DC_V2", "GND_TO_DC_V3", "AC_V12", "AC_V23"]), ...
        get_ASR_samples(1,1,1, 1,1,1, 1,1,1))

      % ==========================
      % BDM = Unknown, DLR Unknown
      % ==========================
      bltsSamplesAVolt(:, :, 1) = AsidTestSamplesSrm(A("DC_V1"));
      bltsSamplesAVolt(:, :, 2) = AsidTestSamplesSrm(A("DC_V2"));
      bltsSamplesAVolt(:, :, 3) = AsidTestSamplesSrm(A("DC_V3"));
      bltsSamplesAVolt(:, :, 4) = AsidTestSamplesSrm(A("AC_V12"));
      bltsSamplesAVolt(:, :, 5) = AsidTestSamplesSrm(A("AC_V23"));
      test(NaN, NaN, ...
        bltsSamplesAVolt, ...
        R([ ...
        "UNKNOWN_TO_NOWHERE", ...
        "UNKNOWN_TO_NOWHERE", ...
        "UNKNOWN_TO_NOWHERE", ...
        "UNKNOWN_TO_NOWHERE", ...
        "AC_V23"]), ...
        get_ASR_samples(0,0,0, 0,0,0, 0,0,1))
    end



    function test_reconstruct_ASR_samples_subsequence(testCase)

      A = bicas.proc.L1L2.const.C.ASID_DICT;



      % Local utility function.
      function assert_relation(A, B, C)
        % NOTE: Uses testCase. ==> Do not make static function.
        b = ~isnan(A) & ~isnan(B) & ~isnan(C);

        testCase.verifyEqual( A(b), B(b) + C(b) )
      end



      % NOTE: Only verifies the correct relationships between the return
      % results. Does not verify entire return results.
      % BUG/NOTE: Will fail if function returns NaN when it should not!
      function test(inputFieldsCa)
        nRows = size(inputFieldsCa{2}, 1);
        AsSrm = bicas.utils.SameRowsMap("uint8", nRows, 'EMPTY');
        for i = 1:(numel(inputFieldsCa)/2)
          asid   = inputFieldsCa{2*i-1};
          sample = inputFieldsCa{2*i  };
          AsSrm.add(asid, sample)
        end

        % RUN FUNCTION TO BE TESTED
        bicas.proc.L1L2.demuxer.reconstruct_ASR_samples_subsequence(AsSrm);
        ActAsSrm = AsSrm;

        % Test all possible relationsships.
        %
        % NOTE: Implicitly asserts that all fields are present.
        % NOTE: Must account for that some fields may be NaN, and
        %       therefore can not be checked against relations.
        assert_relation(ActAsSrm(A("DC_V1")), ActAsSrm(A("DC_V12")), ActAsSrm(A("DC_V2")))
        assert_relation(ActAsSrm(A("DC_V2")), ActAsSrm(A("DC_V23")), ActAsSrm(A("DC_V3")))
        assert_relation(ActAsSrm(A("DC_V1")), ActAsSrm(A("DC_V13")), ActAsSrm(A("DC_V3")))

        % DC. All diffs
        assert_relation(ActAsSrm(A("DC_V13")), ActAsSrm(A("DC_V12")), ActAsSrm(A("DC_V23")))

        % AC. All diffs
        assert_relation(ActAsSrm(A("AC_V13")), ActAsSrm(A("AC_V12")), ActAsSrm(A("AC_V23")))
      end



      test({A("DC_V1"), 19, A("DC_V12"), 27, A("DC_V23"), 33,    A("AC_V12"), 54, A("AC_V23"), 75});    % bdm=0, dlr=0
      test({A("DC_V1"), 19, A("DC_V13"), 27, A("DC_V23"), 33,    A("AC_V13"), 54, A("AC_V23"), 75});    % bdm=0, dlr=1
      test({A("DC_V2"), 19, A("DC_V3"),  27, A("DC_V23"), 19-27, A("AC_V12"), 54, A("AC_V23"), 75});    % bdm=1
      test({A("DC_V1"),  2, A("DC_V2"),   7, A("DC_V3"),  32,    A("AC_V12"), 74, A("AC_V23"), 85});    % bdm=4

    end



    function test_reconstruct_missing_data(testCase)

      % Test data with relationship a1 + a2 == a3 for non-NaN.
      function test_sum(ACa, expACa)
        [actA1,actA2,actA3] = bicas.proc.L1L2.demuxer.reconstruct_missing_data( ...
          ACa{1},        ACa{2},        ACa{3}, ...
          isnan(ACa{1}), isnan(ACa{2}), isnan(ACa{3}), ...
          @(a2, a3) (a3-a2), ...
          @(a1, a3) (a3-a1), ...
          @(a1, a2) (a1+a2));

        testCase.assertEqual(actA1, expACa{1})
        testCase.assertEqual(actA2, expACa{2})
        testCase.assertEqual(actA3, expACa{3})
      end

      %==============
      % Empty arrays
      %==============
      for sizeCa = {[0,0], [1,0], [0,1], [0,0,1]}
        A = zeros(sizeCa{1});
        test_sum({A, A, A}, {A, A, A})
      end

      %====================
      % Scalars, all cases
      %====================
      % Inconsistent existing relationship.
      test_sum({1, 3, 9}, {1, 3, 9})

      % Reconstruct
      test_sum({nan, 3, 5}, {2, 3, 5})
      test_sum({2, nan, 5}, {2, 3, 5})
      test_sum({2, 3, nan}, {2, 3, 5})

      % Can not reconstruct values.
      test_sum({2, nan, nan}, {2, nan, nan})
      test_sum({nan, 3, nan}, {nan, 3, nan})
      test_sum({nan, nan, 5}, {nan, nan, 5})
      %
      test_sum({nan, nan, nan}, {nan, nan, nan})

      %==========================================
      % Non-scalar array, multiple cases at once
      %==========================================
      test_sum({ ...
        [1, nan,   4;   7, nan, nan], ...
        [3,   3, nan;   8,   9, nan], ...
        [9,   5,  10; nan, nan, nan] ...
        }, { ...
        [1,   2,   4;   7, nan, nan], ...
        [3,   3,   6;   8,   9, nan], ...
        [9,   5,  10;  15, nan, nan] ...
        })
    end



    function test_reconstruct_ASR_samples2(testCase)
      N = NaN;

      SAMPLES_AR_DATA = [...
        1,2,4,   -1,-3,-2,   7,10,3;
        1,N,N,   -1, N,-2,   7,10,N;
        1,N,N,   -1, N,-2,   7, N,N;
        ];
      VSTB_AR_DATA = logical([ ...
        1,0,0,    0, 0, 0,   1, 0,0;
        1,0,0,    0, 0, 0,   1, 0,0;
        0,0,0,    0, 0, 1,   1, 0,0;
        ]);
      SdcdDict = testCase.create_SdcdDict(SAMPLES_AR_DATA, VSTB_AR_DATA);



      EXP_SAMPLES_AR_DATA = [...
        1,2,4,   -1,-3,-2,   7,10,3;
        1,2,4,   -1,-3,-2,   7,10,3;
        1,2,4,   -1,-3,-2,   7, N,N;
        ];
      EXP_VSTB_AR_DATA = logical([ ...
        1,0,0,    0, 0, 0,   1, 0,0;
        1,1,1,    0, 0, 0,   1, 0,1;
        0,0,1,    0, 1, 1,   1, 0,0;
        ]);
      ExpSdcdDict = testCase.create_SdcdDict(EXP_SAMPLES_AR_DATA, EXP_VSTB_AR_DATA);



      ActSdcdDict = bicas.proc.L1L2.demuxer.reconstruct_ASR_samples2(SdcdDict);

      % IMPLEMENTATION NOTE: Not comparing entire
      % bicas.proc.L1L2.SdChannelDataDict objects since
      % (1) it will fail also when the objects are identical (since has not
      %     implemented support for testing equality?), and
      % (2) it helps to compare object components separately when debugging
      %     tests.
      for sdid = bicas.proc.L1L2.const.C.SDID_ASR_AR'
        ActSdcd = ActSdcdDict.get(sdid);
        ExpSdcd = ExpSdcdDict.get(sdid);

        if ~isequaln(ActSdcd.samplesAr, ExpSdcd.samplesAr)
          sdid
          ActSdcd.samplesAr
          ExpSdcd.samplesAr
        end
        if ~isequaln(ActSdcd.vstbAr, ExpSdcd.vstbAr)
          sdid
          ActSdcd.vstbAr
          ExpSdcd.vstbAr
        end

        % Check everything (partially overlapping with above).
        testCase.assertEqual(ActSdcd, ExpSdcd)
      end
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Fast-and-easy function for creating one bicas.proc.L1L2.SdChannelDataDict
    % from variables on a format suitable for hardcoding (CWF only).
    function SdcdDict = create_SdcdDict(samplesArData, vstbArData)
      SDID_AR = bicas.proc.L1L2.const.C.SDID_ASR_AR;

      SdcdDict = bicas.proc.L1L2.SdChannelDataDict();
      for i = 1:numel(SDID_AR)
        Sdcd = bicas.proc.L1L2.SdChannelData(...
          samplesArData(:, i), ...
          vstbArData(   :, i));
        SdcdDict = SdcdDict.set(SDID_AR(i), Sdcd);
      end
    end



  end



end
