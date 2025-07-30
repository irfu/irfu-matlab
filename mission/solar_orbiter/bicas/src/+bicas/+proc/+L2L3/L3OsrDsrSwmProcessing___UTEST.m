%
% matlab.unittest automatic test code for bicas.proc.L2L3.L3OsrDsrSwmProcessing.
%
% NOTE: Very low code coverage.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08, from older test code.
%
classdef L3OsrDsrSwmProcessing___UTEST < matlab.unittest.TestCase
  % PROBLEM:
  %   Tests of process_L2_to_L3() do not test QUALITY_FLAG lowered by NSO
  %   table. This results in output QUALITY_FLAG=max, also for science data
  %   blanked by L2-->L3 processing. Input QUALITY_FLAG values are set low but
  %   independently of NSO table and L2QBM which means they can not be
  %   "reconstructed" from NSO table, which means QUALITY_FLAG is re-derived as
  %   high.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_ASSUMPTION(testCase)
      % Tests are designed for this value.
      testCase.fatalAssertTrue(bicas.const.N_MIN_OSR_SAMPLES_PER_BIN == 3)
    end



    % Conceivable special cases for bins to test, including combinations
    % thereof.
    function test_process_L2_to_L3___0(testCase)
      % NOTE: Unclear how much testing is meaningful. Could add more tests.
      %
      % PROPOSITION: The complexity of the test code implies that the
      %              underlying code needs to be refactored somehow.
      %   NOTE: Test sets GAs!
      %   PRO: This test code really tests
      %       bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template() to a large extent.
      %       That function sets:
      %         Epoch
      %         QUALITY_FLAG
      %         QUALITY_BITMASK
      %         L2_QUALITY_BITMASK
      %         DELTA_PLUS_MINUS
      %   PROPOSAL: Reorg into test of
      %       bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template().
      %
      % PROPOSAL: Separate tests (function calls) for different special
      %           cases of bins.
      %   PRO: Easier to follow behaviour in tested code.
      %
      % PROPOSAL: Check OSR data for NaN (not just DSR).

      % Defines bins representing different special cases.
      %
      % NOTE: No consistent relationship between VDC and EDC values, since
      % that is not needed for testing.
      DATA_OSR = [...
        % Too few records in one bin.
        10,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        % Zero records in one bin.
        % - (no data/rows)
        % OBSOLETE SCENARIO: One QUALITY_FLAG is too low, but still enough records for one bin.
        30,   1,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        31,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        32,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        33,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        % OBSOLETE SCENARIO: One QUALITY_FLAG is too low, and therefore NOT enough records for one bin.
        40,   1,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        41,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        42,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        % OBSOLETE SCENARIO: All QUALITY_FLAG are too low for bin.
        50,   1,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        51,   1,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        52,   1,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        ];

      % NOTE: nanData is almost the same as QUALITY_FLAG < 2, except for bin
      % filled with only QUALITY_FLAG=FV (sic!). Might change implementation
      % w.r.t. this behaviour some day.
      DATA_DSR = [...
        10,    4,     0,     0, 1; ...
        20,  255, 65535, 65535, 1; ...
        30,    4,     0,     0, 0; ...
        40,    4,     0,     0, 0; ...
        50,    4,     0,     0, 0; ...
        ]
      bicas.proc.L2L3.L3OsrDsrSwmProcessing___UTEST.test(testCase, ...
        osrIn_Epoch_sec          =DATA_OSR(:, 1), ...
        osrIn_QUALITY_FLAG       =DATA_OSR(:, 2), ...
        osrExp_QUALITY_FLAG      =DATA_OSR(:, 3), ...
        osrIn_QUALITY_BITMASK    =DATA_OSR(:, 4), ...
        osrIn_L2_QUALITY_BITMASK =DATA_OSR(:, 5), ...
        osrIn_VDC                =DATA_OSR(:, 6:8), ...
        osrIn_EDC                =DATA_OSR(:, 9:11), ...
        ...
        dsrExp_Epoch_sec         =DATA_DSR(:, 1), ...
        dsrExp_QUALITY_FLAG      =DATA_DSR(:, 2), ...
        dsrExp_QUALITY_BITMASK   =DATA_DSR(:, 3), ...
        dsrExp_L2_QUALITY_BITMASK=DATA_DSR(:, 4), ...
        dsrExp_nanData           =DATA_DSR(:, 5))
    end



    % Conceivable special cases for bins to test, including combinations
    % thereof.
    function test_process_L2_to_L3___1(testCase)
      N = NaN;

      DATA_OSR = [...
        % One QUALITY_FLAG==FV, but there are still enough QUALITY_FLAG<>FV
        % records for one bin.
        60, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
        61,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        62,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        63,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        % One QUALITY_FLAG==FV, NOT enough other records, enough DATA.
        % ==> Ambiguous.
        % Should never have input QUALITY_FLAG==FV + non-NaN data in
        % the first place.
        70, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
        71,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        72,   2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        % All QUALITY_FLAG==FV
        80, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
        81, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
        82, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
        % All data NaN. ==> Not enough data for one bin.
        90,   2,   4, 0,  0,   N,N,N,  N, N, N; ...
        91,   2,   4, 0,  0,   N,N,N,  N, N, N; ...
        92,   2,   4, 0,  0,   N,N,N,  N, N, N; ...
        % Some data NaN. ==> Not enough records/data for one bin.
        100,  2,   4, 0,  0,   N,N,N,  N, N, N; ...
        101,  2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        102,  2,   4, 0,  0,   1,2,3, -1,-2,-3; ...
        ];
      DATA_DSR = [...
        60,    4,     0,     0, 0; ...
        70,    4,     0,     0, 0; ...
        80,  255,     0,     0, 0; ...
        90,    4,     0,     0, 1; ...
        100,   4,     0,     0, 1; ...
        ];
      bicas.proc.L2L3.L3OsrDsrSwmProcessing___UTEST.test(testCase, ...
        osrIn_Epoch_sec          =DATA_OSR(:, 1), ...
        osrIn_QUALITY_FLAG       =DATA_OSR(:, 2), ...
        osrExp_QUALITY_FLAG      =DATA_OSR(:, 3), ...
        osrIn_QUALITY_BITMASK    =DATA_OSR(:, 4), ...
        osrIn_L2_QUALITY_BITMASK =DATA_OSR(:, 5), ...
        osrIn_VDC                =DATA_OSR(:, 6:8), ...
        osrIn_EDC                =DATA_OSR(:, 9:11), ...
        ...
        dsrExp_Epoch_sec         =DATA_DSR(:, 1), ...
        dsrExp_QUALITY_FLAG      =DATA_DSR(:, 2), ...
        dsrExp_QUALITY_BITMASK   =DATA_DSR(:, 3), ...
        dsrExp_L2_QUALITY_BITMASK=DATA_DSR(:, 4), ...
        dsrExp_nanData           =DATA_DSR(:, 5))
    end



    % Normal bin
    % Test merging QUALITY_BITMASK bits.
    function test_process_L2_to_L3___2(testCase)
      assert(bicas.const.Q.L2_CHANNEL_SATURATION_QRCSM.get("SATURATION_ZV_V12").L2_QUALITY_BITMASK == uint16( 8));
      assert(bicas.const.Q.L2_CHANNEL_SATURATION_QRCSM.get("SATURATION_ZV_V13").L2_QUALITY_BITMASK == uint16(16));
      assert(bicas.const.Q.L2_CHANNEL_SATURATION_QRCSM.get("SATURATION_ZV_V23").L2_QUALITY_BITMASK == uint16(32));

      DATA_OSR = [ ...
        -1,  2, 4,  1,  8,   1,2,3, -1,-2,-3; ...
        0,   2, 4,  2, 16,   1,2,3, -1,-2,-3; ...
        1,   2, 4,  4, 32,   1,2,3, -1,-2,-3; ...
        ];
      % CHANNEL_SATURATION
      DATA_DSR = [ ...
        0,   4,     7,    32+16+8, 0; ...
        ];

      bicas.proc.L2L3.L3OsrDsrSwmProcessing___UTEST.test(testCase, ...
        osrIn_Epoch_sec          =DATA_OSR(:, 1), ...
        osrIn_QUALITY_FLAG       =DATA_OSR(:, 2), ...
        osrExp_QUALITY_FLAG      =DATA_OSR(:, 3), ...
        osrIn_QUALITY_BITMASK    =DATA_OSR(:, 4), ...
        osrIn_L2_QUALITY_BITMASK =DATA_OSR(:, 5), ...
        osrIn_VDC                =DATA_OSR(:, 6:8), ...
        osrIn_EDC                =DATA_OSR(:, 9:11), ...
        ...
        dsrExp_Epoch_sec         =DATA_DSR(:, 1), ...
        dsrExp_QUALITY_FLAG      =DATA_DSR(:, 2), ...
        dsrExp_QUALITY_BITMASK   =DATA_DSR(:, 3), ...
        dsrExp_L2_QUALITY_BITMASK=DATA_DSR(:, 4), ...
        dsrExp_nanData           =DATA_DSR(:, 5))
    end



  end    % methods(Test)



  methods(Static)


    % NOTE: Does not test (L3 DENSITY) L3_QUALITY_BITMASK.
    % NOTE: Does not test science data except for blanking.
    % NOTE: Relies on solo.vdccal() and solo.psp2ne() and that they process
    %       NaN-->NaN.
    %
    function test(testCase, A)
      % PROPOSAL: Have tests explicitly specify the EXCD output and use
      %           bicas.proc.L2L3.ExternalCodeTest.
      %             DCE_SRF_out
      %             PSP_out
      %             ScPot_out
      %             NeScp
      %             NeScpQualityBit
      %   PRO: More rigorous.
      %     PRO: If, not then relying on that solo.vdccal() and solo.psp2ne()
      %   CON: Too much hardcoded input data.
      %     PRO: 3+4=7 more OSR columns.
      %   CON: Makes the submitted EXCD input data redundant.
      %     CON-PROPOSAL: bicas.proc.L2L3.ExternalCodeTest assert the EXCD
      %                   input.
      %       NOTE: EXCD input is TSeries objects for which equality is
      %             uncertain.
      %
      % NOTE: Ideally/theoretically, a test should separately specify
      %   (1) all input+output for
      %       [...
      %         OutEfieldOsr,  OutEfieldDsr, ...
      %         OutScpotOsr,   OutScpotDsr, ...
      %         OutDensityOsr, OutDensityDsr ...
      %       ] = bicas.proc.L2L3.L3OsrDsrSwmProcessing.process_L2_to_L3(...
      %             InLfrCwf, NsoTable, Excd, Bso, L);
      %       ==> Input: 9 columns OSR
      %           Output 2 columns OSR, 4+14=18 columns DSR
      %   (2) all input (for asserting argument values) and output for EXCD
      %       (bicas.proc.L2L3.ExternalCodeTest).
      %       [EdcSrfTs, PspTs, ScpotTs, vdccalCodeVerStr, vdccalMatVerStr] = Excd.vdccal(VdcTs, EdcTs, []);
      %       [NeScpTs, NeScpQualityBitTs, psp2neCodeVerStr]                = Excd.psp2ne(PspTs);
      %       ==> Input:  3+3+1=7 columns OSR
      %           Output: 3+2  =5 columns OSR
      arguments
        testCase
        A.osrIn_Epoch_sec
        A.osrIn_QUALITY_FLAG
        A.osrExp_QUALITY_FLAG
        A.osrIn_QUALITY_BITMASK
        A.osrIn_L2_QUALITY_BITMASK
        A.osrIn_VDC
        A.osrIn_EDC
        %
        A.dsrExp_Epoch_sec
        A.dsrExp_QUALITY_FLAG
        A.dsrExp_QUALITY_BITMASK
        A.dsrExp_L2_QUALITY_BITMASK
        A.dsrExp_nanData
      end

      L   = bicas.Logger('NO_STDOUT', false);
      Bso = bicas.create_default_BSO();
      Bso.override_value('PROCESSING.ZV_QUALITY_FLAG_MAX', 2, mfilename)
      Bso.make_read_only();

      % Tests are designed for this value.
      assert(Bso.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX') == 2)

      FV_QUALITY_FLAG       = uint8(255);
      FV_QUALITY_BITMASK    = uint16(65535);
      FV_L2_QUALITY_BITMASK = uint16(65535);

      % PROPOSAL: Add NSO events.
      NsoTable = bicas.NsoTable(int64.empty(0, 1), int64.empty(0, 1), string.empty(0, 1));

      BASE_TT2000 = spdfparsett2000('2020-03-14T00:00:00');

      % Input OSR
      InLfrCwf.Ga.OBS_ID    = {' '};
      InLfrCwf.Ga.SOOP_TYPE = {' '};
      InLfrCwf.Zv.Epoch                 = int64(                     A.osrIn_Epoch_sec*1e9) + BASE_TT2000;
      InLfrCwf.ZvFpa.QUALITY_FLAG       = bicas.utils.FPArray(uint8( A.osrIn_QUALITY_FLAG),       'FILL_VALUE', FV_QUALITY_FLAG);
      InLfrCwf.ZvFpa.QUALITY_BITMASK    = bicas.utils.FPArray(uint16(A.osrIn_QUALITY_BITMASK),    'FILL_VALUE', FV_QUALITY_BITMASK);
      InLfrCwf.ZvFpa.L2_QUALITY_BITMASK = bicas.utils.FPArray(uint16(A.osrIn_L2_QUALITY_BITMASK), 'FILL_VALUE', FV_L2_QUALITY_BITMASK);
      InLfrCwf.ZvFpa.DELTA_PLUS_MINUS   = bicas.utils.FPArray(int64(ones(size(InLfrCwf.Zv.Epoch))) * mode(diff(InLfrCwf.Zv.Epoch)));
      InLfrCwf.ZvFpa.VDC                = bicas.utils.FPArray(A.osrIn_VDC, 'FILL_VALUE', NaN).cast('single');
      InLfrCwf.ZvFpa.EDC                = bicas.utils.FPArray(A.osrIn_EDC, 'FILL_VALUE', NaN).cast('single');

      % Expected OSR
      ExpOsr.Zv.QUALITY_FLAG       = bicas.utils.FPArray(uint8(A.osrExp_QUALITY_FLAG),        'FILL_VALUE', FV_QUALITY_FLAG);
      % Expected DSR
      ExpDsr.Zv.Epoch              = int64(                     A.dsrExp_Epoch_sec*1e9) + BASE_TT2000;
      ExpDsr.Zv.QUALITY_FLAG       = bicas.utils.FPArray(uint8( A.dsrExp_QUALITY_FLAG),       'FILL_VALUE', FV_QUALITY_FLAG);
      ExpDsr.Zv.QUALITY_BITMASK    = bicas.utils.FPArray(uint16(A.dsrExp_QUALITY_BITMASK),    'FILL_VALUE', FV_QUALITY_BITMASK);
      ExpDsr.Zv.L2_QUALITY_BITMASK = bicas.utils.FPArray(uint16(A.dsrExp_L2_QUALITY_BITMASK), 'FILL_VALUE', FV_L2_QUALITY_BITMASK);
      ExpDsr.nanData               = logical(A.dsrExp_nanData);

      % ================
      % Set RVs for EXCD
      % ================
      % IMPLEMENTATION NOTE: Currently setting all the science data to nonsense
      % and not checking the output result. Pure downsampling should be tested
      % elsewhere.
      % vectorAr = ones(size(InLfrCwf.ZvFpa.VDC      ));   % Time series of vectors.
      % vectorAr(:, 1) = 0;                                % E-field must always have zero X component.
      % scalarAr = ones(size(InLfrCwf.ZvFpa.VDC(:, 1)));   % Time series of scalars.
      %
      % VdccalRv = [];
      % VdccalRv.DCE_SRF_out                  = TSeries(...
      %   EpochTT(InLfrCwf.Zv.Epoch), ...
      %   vectorAr, ...
      %   'TensorOrder', 1, ...
      %   'repres',      {'x', 'y', 'z'});
      % VdccalRv.DCE_SRF_out.units            = 'mV/m';
      % VdccalRv.DCE_SRF_out.coordinateSystem = 'SRF';
      % VdccalRv.PSP_out                      = TSeries(...
      %   EpochTT(InLfrCwf.Zv.Epoch), ...
      %   scalarAr, ...
      %   'TensorOrder', 0);
      % VdccalRv.PSP_out.units                = 'V';
      % VdccalRv.ScPot_out                    = TSeries(...
      %   EpochTT(InLfrCwf.Zv.Epoch), ...
      %   scalarAr, ...
      %   'TensorOrder', 0);
      % VdccalRv.ScPot_out.units              = 'V';
      % VdccalRv.codeVerStr                   = '2022-12-06T13:23:14';
      % VdccalRv.matVerStr                    = 'd23K123_20230707.mat';
      % %-------------------------------------------------------------------
      % Psp2neRv = [];
      % Psp2neRv.NeScp           = TSeries(...
      %   EpochTT(InLfrCwf.Zv.Epoch), scalarAr, ...
      %   'TensorOrder', 0);
      % Psp2neRv.NeScp.units     = 'cm^-3';
      % Psp2neRv.NeScpQualityBit = TSeries(...
      %   EpochTT(InLfrCwf.Zv.Epoch), scalarAr, ...
      %   'TensorOrder', 0);
      % Psp2neRv.codeVerStr      = '2023-08-11T10:11:00';

      Excd = bicas.proc.L2L3.ExternalCodeImplementation();
      % Excd = bicas.proc.L2L3.ExternalCodeTest(VdccalRv, Psp2neRv);
      % IMPLEMENTATION NOTE: Does not yet use bicas.proc.L2L3.ExternalCodeTest
      % since it is quite complicated to specify the output.
      % NOTE: Not using bicas.proc.L2L3.ExternalCodeTest means relying on that
      % solo.vdccal() and solo.psp2ne() always convert NaN-->NaN for the same
      % record.

      %##################################################################
      % CALL CODE TO BE TESTED
      %##################################################################
      [ActEfieldOsr,  ActEfieldDsr, ...
        ActScpotOsr,   ActScpotDsr, ...
        ActDensityOsr, ActDensityDsr] ...
        = bicas.proc.L2L3.L3OsrDsrSwmProcessing.process_L2_to_L3(InLfrCwf, NsoTable, Excd, Bso, L);
      %##################################################################

      % OSR
      for ActOsrCa = {ActEfieldOsr, ActScpotOsr, ActDensityOsr}'
        ActOsr = ActOsrCa{1}.Zv;
        testCase.assertEqual(ActOsr.Epoch,              InLfrCwf.Zv.Epoch)
        testCase.assertEqual(ActOsr.QUALITY_FLAG,       ExpOsr.Zv.QUALITY_FLAG)
        testCase.assertEqual(ActOsr.QUALITY_BITMASK,    InLfrCwf.ZvFpa.QUALITY_BITMASK)
        testCase.assertEqual(ActOsr.L2_QUALITY_BITMASK, InLfrCwf.ZvFpa.L2_QUALITY_BITMASK)
      end

      % DSR
      for ActDsrCa = {ActEfieldDsr, ActScpotDsr, ActDensityDsr}'
        ActDsr = ActDsrCa{1}.Zv;
        testCase.assertEqual(ActDsr.Epoch,              ExpDsr.Zv.Epoch)
        testCase.assertEqual(ActDsr.QUALITY_FLAG,       ExpDsr.Zv.QUALITY_FLAG)
        testCase.assertEqual(ActDsr.QUALITY_BITMASK,    ExpDsr.Zv.QUALITY_BITMASK)
        testCase.assertEqual(ActDsr.L2_QUALITY_BITMASK, ExpDsr.Zv.L2_QUALITY_BITMASK)
      end

      testCase.assertEqual(all(ActEfieldDsr.Zv.EDC_SRF.fpAr,    2), ExpDsr.nanData);
      testCase.assertEqual(all(ActEfieldDsr.Zv.EDCSTD_SRF.fpAr, 2), ExpDsr.nanData);
      testCase.assertEqual(    ActDensityDsr.Zv.DENSITY.fpAr      , ExpDsr.nanData);
      testCase.assertEqual(    ActDensityDsr.Zv.DENSITYSTD.fpAr   , ExpDsr.nanData);
      testCase.assertEqual(    ActScpotDsr.Zv.SCPOT.fpAr          , ExpDsr.nanData);
      testCase.assertEqual(    ActScpotDsr.Zv.SCPOTSTD.fpAr       , ExpDsr.nanData);
      testCase.assertEqual(    ActScpotDsr.Zv.PSP.fpAr            , ExpDsr.nanData);
      testCase.assertEqual(    ActScpotDsr.Zv.PSPSTD.fpAr         , ExpDsr.nanData);
    end



  end    % methods(Static)



end
