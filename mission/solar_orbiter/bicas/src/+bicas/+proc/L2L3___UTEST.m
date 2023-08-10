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
            % QUALITY_FLAG = fill value
            % #samples < bicas.constants.N_MIN_SAMPLES_PER_DSR_BIN
            % 
            % NOTE: Makes no sense testing the actual output data, since it is
            % dependent on BICAS-external functions, except maybe whether it is
            % NaN or not.

            % Test(s) are designed for this value.
            assert(bicas.constants.N_MIN_SAMPLES_PER_DSR_BIN == 3)



            N = NaN;
            L        = bicas.Logger('none', false);
            SETTINGS = bicas.create_default_SETTINGS();
            SETTINGS.make_read_only();



            InLfrCwf.ZvFv = struct(...
                'QUALITY_FLAG',       uint8(255), ...
                'QUALITY_BITMASK',    uint16(65535), ...
                'L2_QUALITY_BITMASK', uint16(65535));
            InLfrCwf.Ga.OBS_ID    = {' '};
            InLfrCwf.Ga.SOOP_TYPE = {' '};

            %===================================================================
            % ORIS DATA: Input+expected output
            % --------------------------------
            % Defines bins representing different special cases.
            %
            % NOTE: No consistent relationship between VDC and EDC values, since
            % that is not needed for testing.
            %===================================================================
            % Columns: Epoch [s], QUALITY_FLAG_in, QUALITY_FLAG_ORIS_out, ...
            DATA1 = [...
                % Normal bin
                % Test merging QUALITY_BITMASK bits.
                 -1,   2,   2, 1,  8,   1,2,3, -1,-2,-3; ...
                  0,   2,   2, 2, 16,   1,2,3, -1,-2,-3; ...
                  1,   2,   2, 4, 32,   1,2,3, -1,-2,-3; ...
                % Too few records
                 10,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                % Zero records
                % - (no data)
                % One QUALITY_FLAG too low, but still enough records.
                 30,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 31,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 32,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 33,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                % One QUALITY_FLAG too low, and therefore NOT enough records.
                 40,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 41,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 42,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                % All QUALITY_FLAG too low.
                 50,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 51,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 52,   1, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 % One QUALITY_FLAG==fill value, but still enough records.
                 60, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 61,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 62,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 63,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 % One QUALITY_FLAG==fill value, NOT enough other records, enough DATA.
                 % ==> Ambiguous.
                 % Should never have input QUALITY_FLAG==fill value + non-NaN data in
                 % the first place.
                 70, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 71,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 72,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                 % All QUALITY_FLAG==fill value
                 80, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 81, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                 82, 255, 255, 0,  0,   1,2,3, -1,-2,-3; ...
                % All data NaN
                 90,   2, 255, 0,  0,   N,N,N, N,N,N; ...
                 91,   2, 255, 0,  0,   N,N,N, N,N,N; ...
                 92,   2, 255, 0,  0,   N,N,N, N,N,N; ...
                % Some data NaN, ==> Not enough records/data
                100,   2, 255, 0,  0,   N,N,N, N,N,N; ...
                101,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                102,   2,   2, 0,  0,   1,2,3, -1,-2,-3; ...
                ];

            %===========
            % DSR DATA
            %===========
            DATA2 = [...
                  0,   2,     7,    56; ...
                 10, 255,     0,     0; ...
                 20, 255, 65535, 65535; ...
                 30,   2,     0,     0; ...
                 40, 255,     0,     0; ...
                 50, 255,     0,     0; ...
                 60,   2,     0,     0; ...
                 70,   2,     0,     0; ...
                 80, 255,     0,     0; ...
                 90, 255,     0,     0; ...
                100, 255,     0,     0; ...
                ];

            BASE_TT2000 = spdfparsett2000('2020-03-14T00:00:00');

            InLfrCwf.Zv.Epoch              = int64( DATA1(:, 1)*1e9) + BASE_TT2000;
            InLfrCwf.Zv.QUALITY_FLAG       = uint8( DATA1(:, 2));
            InLfrCwf.Zv.QUALITY_BITMASK    = uint16(DATA1(:, 4));
            InLfrCwf.Zv.L2_QUALITY_BITMASK = uint16(DATA1(:, 5));
            InLfrCwf.Zv.DELTA_PLUS_MINUS   = int64(ones(size(InLfrCwf.Zv.Epoch))) * mode(diff(InLfrCwf.Zv.Epoch));
            InLfrCwf.Zv.VDC = single(DATA1(:, 6: 8));
            InLfrCwf.Zv.EDC = single(DATA1(:, 9:11));

            ExpOris.Zv.QUALITY_FLAG       = uint8( DATA1(:, 3));
            %
            ExpDsr.Zv.Epoch              = int64( DATA2(:, 1)*1e9) + BASE_TT2000;
            ExpDsr.Zv.QUALITY_FLAG       = uint8( DATA2(:, 2));
            ExpDsr.Zv.QUALITY_BITMASK    = uint16(DATA2(:, 3));
            ExpDsr.Zv.L2_QUALITY_BITMASK = uint16(DATA2(:, 4));

            [OutEfieldOris,  OutEfieldDsr, ...
             OutScpotOris,   OutScpotDsr, ...
             OutDensityOris, OutDensityDsr] ...
            = bicas.proc.L2L3.process_L2_to_L3(InLfrCwf, SETTINGS, L);

            % ORIS
            testCase.verifyEqual(OutEfieldOris.Zv.Epoch,              InLfrCwf.Zv.Epoch)
            testCase.verifyEqual(OutEfieldOris.Zv.QUALITY_FLAG,       ExpOris.Zv.QUALITY_FLAG)
            testCase.verifyEqual(OutEfieldOris.Zv.QUALITY_BITMASK,    InLfrCwf.Zv.QUALITY_BITMASK)
            testCase.verifyEqual(OutEfieldOris.Zv.L2_QUALITY_BITMASK, InLfrCwf.Zv.L2_QUALITY_BITMASK)

            % DSR
            testCase.verifyEqual(OutEfieldDsr.Zv.Epoch,              ExpDsr.Zv.Epoch)
            testCase.verifyEqual(OutEfieldDsr.Zv.QUALITY_FLAG,       ExpDsr.Zv.QUALITY_FLAG)
            testCase.verifyEqual(OutEfieldDsr.Zv.QUALITY_BITMASK,    ExpDsr.Zv.QUALITY_BITMASK)
            testCase.verifyEqual(OutEfieldDsr.Zv.L2_QUALITY_BITMASK, ExpDsr.Zv.L2_QUALITY_BITMASK)
        end
        
        
        
    end    % methods(Test)

    
    
end
