%
% matlab.unittest automatic test code for bicas.proc.L2L3.
%
% NOTE: Very low code coverage.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08, from older test code.
%
classdef L2L3___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        function test_process_L2_to_L3(testCase)
            % Conceivable special cases for bins to test, including combinations
            % thereof
            % ------------------------------------------------------------------
            % No records     (can not be combined with other cases)
            % Only NaN data
            % Partially NaN data
            % QUALITY_FLAG < min value
            % QUALITY_FLAG = FV
            % #samples < bicas.const.N_MIN_OSR_SAMPLES_PER_BIN
            % 
            % Unclear how much testing is meaningful. Could add more tests.
            %
            % PROPOSITION: The complexity of the test code implies that the
            %              underlying code needs to be refactored somehow.
            %   NOTE: Test sets GAs!
            %   PRO: This test code really tests
            %       bicas.proc.dsr.init_shared_DSR_ZVs() to a large extent. That
            %       function sets:
            %         Epoch
            %         QUALITY_FLAG
            %         QUALITY_BITMASK
            %         L2_QUALITY_BITMASK
            %         DELTA_PLUS_MINUS
            %   PROPOSAL: Reorg into test of
            %       bicas.proc.dsr.init_shared_DSR_ZVs().
            %
            % PROPOSAL: Separate tests (function calls) for different special
            %           cases of bins.
            %   PRO: Easier to follow behaviour in tested code.
            %
            % PROPOSAL: Check OSR data for NaN (not just DSR).

            % Test(s) are designed for this value.
            assert(bicas.const.N_MIN_OSR_SAMPLES_PER_BIN == 3)



            N = NaN;
            %===================================================================
            % OSR DATA: Input+expected output
            % --------------------------------
            % Defines bins representing different special cases.
            %
            % NOTE: No consistent relationship between VDC and EDC values, since
            % that is not needed for testing.
            %===================================================================
            % Columns:
            % 1:    Epoch [s]
            % 2:    QUALITY_FLAG_in
            % 3:    QUALITY_FLAG_OSR_out
            % 4:    QUALITY_BITMASK_in
            % 5:    L2_QUALITY_BITMASK
            % 6- 8: VDC
            % 9-11: EDC
            DATA_OSR = [...
                % Too few records
                 10,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                % Zero records
                % - (no data)
                % One QUALITY_FLAG is too low, but still enough records.
                 30,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 31,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 32,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 33,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                % One QUALITY_FLAG is too low, and therefore NOT enough records.
                 40,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 41,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 42,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                % All QUALITY_FLAG are too low.
                 50,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 51,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 52,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 % One QUALITY_FLAG==FV, but still enough records.
                 60, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 61,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 62,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 63,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 % One QUALITY_FLAG==FV, NOT enough other records, enough DATA.
                 % ==> Ambiguous.
                 % Should never have input QUALITY_FLAG==FV + non-NaN data in
                 % the first place.
                 70, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 71,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 72,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 % All QUALITY_FLAG==FV
                 80, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 81, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 82, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                % All data NaN
                 90,   2, 255, 0,  0,   N,N,N,  N, N, N; ...
                 91,   2, 255, 0,  0,   N,N,N,  N, N, N; ...
                 92,   2, 255, 0,  0,   N,N,N,  N, N, N; ...
                % Some data NaN, ==> Not enough records/data
                100,   2, 255, 0,  0,   N,N,N,  N, N, N; ...
                101,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                102,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
            ];

            %===========
            % DSR DATA
            %===========
            % Columns:
            % 1: Epoch [s]
            % 2: QUALITY_FLAG
            % 3: QUALITY_BITMASK
            % 4: L2_QUALITY_BITMASK
            % 5: nanData
            % --
            % NOTE: nanData is almost same as QUALITY_FLAG < 2, except for bin
            % filled with only QUALITY_FLAG=FV (sic!). Might change
            % implementation w.r.t. this behaviour some day.
            DATA_DSR = [...
                 10, 255,     0,     0, 1;...
                 20, 255, 65535, 65535, 1; ...
                 30,   2,     0,     0, 0; ...
                 40, 255,     0,     0, 1; ...
                 50, 255,     0,     0, 1; ...
                 60,   2,     0,     0, 0; ...
                 70,   2,     0,     0, 0; ...
                 80, 255,     0,     0, 0; ...
                 90, 255,     0,     0, 1; ...
                100, 255,     0,     0, 1; ...
            ];
            bicas.proc.L2L3___UTEST.test(testCase, DATA_OSR, DATA_DSR)
            
            % Normal bin
            % Test merging QUALITY_BITMASK bits.
            bicas.proc.L2L3___UTEST.test(testCase, ...
                [ ...
                 -1,   2,   2, 1,  8,   1,2,3, -1,-2,-3; ...
                  0,   2,   2, 2, 16,   1,2,3, -1,-2,-3; ...
                  1,   2,   2, 4, 32,   1,2,3, -1,-2,-3; ...
                ], [ ...
                  0,   2,     7,    56, 0; ...
                ] ...
            )

        end
        
        
        
    end    % methods(Test)



    methods(Static)

        % DATA_OSR
        % --------
        % Columns:
        % 1:    Epoch [s]
        % 2:    QUALITY_FLAG_in
        % 3:    QUALITY_FLAG_OSR_out
        % 4:    QUALITY_BITMASK_in
        % 5:    L2_QUALITY_BITMASK
        % 6- 8: VDC
        % 9-11: EDC
        %
        % DATA_DSR
        % --------
        % Columns:
        % 1: Epoch [s]
        % 2: QUALITY_FLAG
        % 3: QUALITY_BITMASK
        % 4: L2_QUALITY_BITMASK
        % 5: nanData
        %
        function test(testCase, DATA_OSR, DATA_DSR)
            L        = bicas.Logger('none', false);
            SETTINGS = bicas.create_default_SETTINGS();
            SETTINGS.make_read_only();

            % Tests are designed for this value.
            assert(SETTINGS.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX') == 2)

            FV_QUALITY_FLAG       = uint8(255);
            FV_QUALITY_BITMASK    = uint16(65535);
            FV_L2_QUALITY_BITMASK = uint16(65535);

            InLfrCwf.Ga.OBS_ID    = {' '};
            InLfrCwf.Ga.SOOP_TYPE = {' '};

            
            
            BASE_TT2000 = spdfparsett2000('2020-03-14T00:00:00');

            InLfrCwf.Zv.Epoch                 = int64( DATA_OSR(:, 1)*1e9) + BASE_TT2000;
            InLfrCwf.ZvFpa.QUALITY_FLAG       = bicas.utils.FPArray(uint8( DATA_OSR(:, 2)), 'FILL_VALUE', FV_QUALITY_FLAG);
            InLfrCwf.ZvFpa.QUALITY_BITMASK    = bicas.utils.FPArray(uint16(DATA_OSR(:, 4)), 'FILL_VALUE', FV_QUALITY_BITMASK);
            InLfrCwf.ZvFpa.L2_QUALITY_BITMASK = bicas.utils.FPArray(uint16(DATA_OSR(:, 5)), 'FILL_VALUE', FV_L2_QUALITY_BITMASK);
            InLfrCwf.ZvFpa.DELTA_PLUS_MINUS   = bicas.utils.FPArray(int64(ones(size(InLfrCwf.Zv.Epoch))) * mode(diff(InLfrCwf.Zv.Epoch)), 'NO_FILL_POSITIONS');
            InLfrCwf.ZvFpa.VDC                = bicas.utils.FPArray(DATA_OSR(:, 6: 8), 'FILL_VALUE', NaN).cast('single');
            InLfrCwf.ZvFpa.EDC                = bicas.utils.FPArray(DATA_OSR(:, 9:11), 'FILL_VALUE', NaN).cast('single');

            ExpOsr.Zv.QUALITY_FLAG       = bicas.utils.FPArray(uint8( DATA_OSR(:, 3)), 'FILL_VALUE', uint8 (FV_QUALITY_FLAG));
            %
            ExpDsr.Zv.Epoch              = int64( DATA_DSR(:, 1)*1e9) + BASE_TT2000;
            ExpDsr.Zv.QUALITY_FLAG       = bicas.utils.FPArray(uint8( DATA_DSR(:, 2)), 'FILL_VALUE', FV_QUALITY_FLAG);
            ExpDsr.Zv.QUALITY_BITMASK    = bicas.utils.FPArray(uint16(DATA_DSR(:, 3)), 'FILL_VALUE', FV_QUALITY_BITMASK);
            ExpDsr.Zv.L2_QUALITY_BITMASK = bicas.utils.FPArray(uint16(DATA_DSR(:, 4)), 'FILL_VALUE', FV_L2_QUALITY_BITMASK);
            ExpDsr.nanData               = logical(DATA_DSR(:, 5));

            % ==========================
            % Set RVs for external code
            % ==========================
            % IMPLEMENTATION NOTE: Currently setting all the science data to
            % nonsense and not checking the output result. Pure downsampling
            % should be tested elsewhere.
            vectorAr = ones(size(InLfrCwf.ZvFpa.VDC      ));   % Time series of vectors.
            scalarAr = ones(size(InLfrCwf.ZvFpa.VDC(:, 1)));   % Time series of scalars.
            %
            VdccalRv = [];
            VdccalRv.DCE_SRF_out = TSeries(...
                EpochTT(InLfrCwf.Zv.Epoch), vectorAr, ...
                'TensorOrder', 1, ...
                'repres',      {'x', 'y', 'z'});
            VdccalRv.PSP_out     = TSeries(...
                EpochTT(InLfrCwf.Zv.Epoch), scalarAr, ...
                'TensorOrder', 0);
            VdccalRv.ScPot_out   = TSeries(...
                EpochTT(InLfrCwf.Zv.Epoch), scalarAr, ...
                'TensorOrder', 0);
            VdccalRv.codeVerStr  = '2022-12-06T13:23:14';
            VdccalRv.matVerStr   = 'd23K123_20230707.mat';
            %-------------------------------------------------------------------
            Psp2neRv = [];
            Psp2neRv.NeScp           = TSeries(...
                EpochTT(InLfrCwf.Zv.Epoch), scalarAr, ...
                'TensorOrder', 0);
            Psp2neRv.NeScpQualityBit = TSeries(...
                EpochTT(InLfrCwf.Zv.Epoch), scalarAr, ...
                'TensorOrder', 0);
            Psp2neRv.codeVerStr      = '2023-08-11T10:11:00';
            
            Ec = bicas.proc.L2L3.ExternalCodeImplementation();
            
            %##################################################################
            % CALL CODE TO BE TESTED
            %##################################################################
            [OutEfieldOsr,  OutEfieldDsr, ...
             OutScpotOsr,   OutScpotDsr, ...
             OutDensityOsr, OutDensityDsr] ...
                = bicas.proc.L2L3.process_L2_to_L3(InLfrCwf, Ec, SETTINGS, L);
            %##################################################################

            % OSR
            for OutOsrCa = {OutEfieldOsr, OutScpotOsr, OutDensityOsr}'
                OutOsr = OutOsrCa{1}.Zv;
                testCase.assertEqual(OutOsr.Epoch,              InLfrCwf.Zv.Epoch)
                testCase.assertEqual(OutOsr.QUALITY_FLAG,       ExpOsr.Zv.QUALITY_FLAG)
                testCase.assertEqual(OutOsr.QUALITY_BITMASK,    InLfrCwf.ZvFpa.QUALITY_BITMASK)
                testCase.assertEqual(OutOsr.L2_QUALITY_BITMASK, InLfrCwf.ZvFpa.L2_QUALITY_BITMASK)
            end

            % DSR
            for OutDsrCa = {OutEfieldDsr, OutScpotDsr, OutDensityDsr}'
                OutDsr = OutDsrCa{1}.Zv;
                testCase.assertEqual(OutDsr.Epoch,              ExpDsr.Zv.Epoch)
                testCase.assertEqual(OutDsr.QUALITY_FLAG,       ExpDsr.Zv.QUALITY_FLAG)
                testCase.assertEqual(OutDsr.QUALITY_BITMASK,    ExpDsr.Zv.QUALITY_BITMASK)
                testCase.assertEqual(OutDsr.L2_QUALITY_BITMASK, ExpDsr.Zv.L2_QUALITY_BITMASK)
            end
            
            testCase.assertEqual(all(OutEfieldDsr.Zv.EDC_SRF.fpAr,    2),   ExpDsr.nanData);
            testCase.assertEqual(all(OutEfieldDsr.Zv.EDCSTD_SRF.fpAr, 2),   ExpDsr.nanData);
            testCase.assertEqual(    OutDensityDsr.Zv.DENSITY.fpAr        , ExpDsr.nanData);
            testCase.assertEqual(    OutDensityDsr.Zv.DENSITYSTD.fpAr     , ExpDsr.nanData);
            testCase.assertEqual(    OutScpotDsr.Zv.SCPOT.fpAr            , ExpDsr.nanData);
            testCase.assertEqual(    OutScpotDsr.Zv.SCPOTSTD.fpAr         , ExpDsr.nanData);
            testCase.assertEqual(    OutScpotDsr.Zv.PSP.fpAr              , ExpDsr.nanData);
            testCase.assertEqual(    OutScpotDsr.Zv.PSPSTD.fpAr           , ExpDsr.nanData);
        end


    
    end    % methods(Static)


    
end
