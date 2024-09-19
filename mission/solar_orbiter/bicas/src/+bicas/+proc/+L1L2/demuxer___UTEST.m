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



  methods(Static, Access=private)



    function AsidTestSamplesSrm = create_channel_test_data(sampleSize)
      % (iChannel, 1) = ASID ID string.
      % (iChannel, 2) = Sample value (that is consistent with other
      %                 channels).
      A = bicas.sconst.C.S_ASID_DICT;
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
        "bicas.proc.L1L2.AntennaSignalId", sampleSize(1), 'EMPTY');
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



    % Test two functions in combination.
    %
    % IMPLEMENTATION NOTE: The design is for historical reasons before the
    % two functions were split up.
    function test_get_routings_calibrated_BLTSs_to_all_ASRs(testCase)

      A = bicas.sconst.C.S_ASID_DICT;
      R = bicas.sconst.C.S_ROUTING_DICT;

      % =========
      % Test data
      % =========
      TEST_DATA_UNKNOWN = [999];   % Data from unknown source.
      SAMPLES_SIZE      = [3,2];

      nRows = SAMPLES_SIZE(1);
      %[channelIdStringsCa, testSamplesCa] = testCase.create_channel_test_data(SAMPLES_SIZE);
      AsidTestSamplesSrm = testCase.create_channel_test_data(SAMPLES_SIZE);



      % Test any BDM.
      function test(bdmFloatNan, dlrFloatNan, bltsSamplesAVolt, ExpRoutingArray, ExpAsrSamplesAVoltSrm)
        assert(numel(ExpRoutingArray) == 5)

        dlrFpa = bicas.utils.FPArray.floatNan2logical(dlrFloatNan);
        bdmFpa = bicas.utils.FPArray.floatNan2int(bdmFloatNan, 'uint8');

        % CALL FUNCTIONS
        ActRoutingArray       = bicas.proc.L1L2.demuxer.get_routings(...
          bdmFpa, dlrFpa);
        ActAsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_all_ASRs(...
          [ActRoutingArray.Sdid], bltsSamplesAVolt);

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
          if Routing.Ssid.is_ASR()
            tempBltsSamplesAVolt(:, :, i) = AsidTestSamplesSrm(Routing.Ssid.Asid);
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

        % Define which varargin{i} corresponds to which ASID.
        ASID_ARRAY = A([
          "DC_V1",  "DC_V2",  "DC_V3",  ...
          "DC_V12", "DC_V13", "DC_V23", ...
          "AC_V12", "AC_V13", "AC_V23" ...
          ]);
        AsSrm = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", nRows, 'EMPTY');

        for iAsid = 1:9
          asidName = ASID_ARRAY(iAsid);

          samplesAVolt = AsidTestSamplesSrm(asidName);
          if ~varargin{iAsid}
            samplesAVolt = samplesAVolt * NaN;
          end
          AsSrm.add(asidName, samplesAVolt);
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

      % ==============
      % BDM = Unknwon, DLR Unknown
      % ==============
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



    function test_complement_ASR(testCase)

      A = bicas.sconst.C.S_ASID_DICT;



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
        AsSrm = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", nRows, 'EMPTY');
        for i = 1:(numel(inputFieldsCa)/2)
          Asid   = inputFieldsCa{2*i-1};
          sample = inputFieldsCa{2*i  };
          AsSrm.add(Asid, sample)
        end

        % RUN FUNCTION TO BE TESTED
        bicas.proc.L1L2.demuxer.complement_ASR(AsSrm);
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



  end    % methods(Test)



end
