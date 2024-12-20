%
% matlab.unittest automatic test code for bicas.proc.L1L2.qual, except for
% method bicas.proc.L1L2.qual.sliding_window_over_fraction().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_set_voltage_current_FV(testCase)

      function test(...
          zv_Epoch, zvUfv, zvAsrSamplesAVoltSrm, zvCurrentAAmpere, ...
          expZvAsrSamplesAVoltSrm, expZvCurrentAAmpere)

        nRows = irf.assert.sizes(...
          zv_Epoch,         [-1], ...
          zvCurrentAAmpere, [-1, 3]);
        assert(zvAsrSamplesAVoltSrm.nRows == nRows)
        L = bicas.Logger('NO_STDOUT', false);

        % NOTE: Modifies argument zvAsrSamplesAVoltSrm (handle object).
        actZvCurrentAAmpere = bicas.proc.L1L2.qual.set_voltage_current_FV(...
          zv_Epoch, zvAsrSamplesAVoltSrm, zvCurrentAAmpere, zvUfv, L);

        actZvAsrSamplesAVoltSrm = zvAsrSamplesAVoltSrm;
        testCase.verifyEqual(actZvAsrSamplesAVoltSrm,  expZvAsrSamplesAVoltSrm)
        testCase.verifyEqual(actZvCurrentAAmpere, expZvCurrentAAmpere)
      end

      %===================================================================

      % Empty data.
      test( ...
        int64(zeros(0, 1)), ...
        false(0, 1), ...
        bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(0, 1)), ...
        zeros(0, 3),  ...
        bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(0, 1)), ...
        zeros(0, 3)  ...
        )

      % Non-empty input data that is not altered.
      test( ...
        int64([10, 11, 12, 13, 14]'), ...
        false(5, 1), ...
        bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(5, 1)), ...
        zeros(5, 3),  ...
        bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(5, 1)), ...
        zeros(5, 3)  ...
        )

      % Non-empty input data that is altered: 1 unaltered + 1 altered
      test( ...
        int64([10, 11]'), ...
        logical([0, 1]'), ...
        bicas.proc.L1L2.qual___UTEST.AsSrm([1, 2]'), ...
        [1:3; 11:13],  ...
        bicas.proc.L1L2.qual___UTEST.AsSrm([1, NaN]'), ...
        [1:3; NaN(1, 3)] ...
        )
      % Non-empty input data that is altered: 2 unaltered + (1+2) altered
      test( ...
        int64([10, 11, 12, 13, 14]'), ...
        logical([0, 1, 0, 1, 1]'), ...
        bicas.proc.L1L2.qual___UTEST.AsSrm([1, 2, 3, 4, 5]'), ...
        [1:3; 11:13; 21:23; 31:33; 41:43],  ...
        bicas.proc.L1L2.qual___UTEST.AsSrm([1, NaN, 3, NaN, NaN]'), ...
        [1:3; NaN(1,3); 21:23; NaN(2,3)] ...
        )

    end



    function test_get_UFV_from_removing_BDMs(testCase)

      function actZvUfv = call_func(zv_Epoch, zvBdmDoubleNan, isLfr, bdmRemoveArray, lfrBdmMarginSec, tdsBdmMarginSec)
        assert(isfloat(zvBdmDoubleNan))
        zvBdmFpa = bicas.utils.FPArray(zvBdmDoubleNan, 'FILL_VALUE', NaN).cast('int8');

        Bso = bicas.create_default_BSO();
        Bso.override_value('PROCESSING.L2.REMOVE_DATA.MUX_MODES',             bdmRemoveArray,  'test')
        Bso.override_value('PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S', lfrBdmMarginSec, 'test')
        Bso.override_value('PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S', tdsBdmMarginSec, 'test')
        Bso.make_read_only()

        L = bicas.Logger('NO_STDOUT', false);

        actZvUfv = bicas.proc.L1L2.qual.get_UFV_from_removing_BDMs(...
          zv_Epoch, zvBdmFpa, isLfr, Bso, L);
      end



      function test(zv_Epoch, zvBdm, isLfr, bdmRemoveArray, lfrBdmMarginSec, tdsBdmMarginSec, expZvUfv)
        irf.assert.sizes(...
          zv_Epoch, [-1, 1], ...
          zvBdm,    [-1, 1], ...
          expZvUfv, [-1, 1] ...
          )

        actZvUfv = call_func(zv_Epoch, zvBdm, isLfr, bdmRemoveArray, lfrBdmMarginSec, tdsBdmMarginSec);

        testCase.verifyEqual(actZvUfv, expZvUfv)
      end
      %===================================================================



      if 1
        for isLfr = [false, true]
          % NOTE: Test PROCESSING.L2.REMOVE_DATA.MUX_MODES = [] (0x0)
          % as is likely to be set to.
          % Iterate over bdmRemoveArray values without BDM=0.
          for bdmRemoveArrayCa = {zeros(1,0), zeros(0,1), [2]', [1,2,3]'}
            bdmRemoveArray = bdmRemoveArrayCa{1};

            % Zero records
            test(...
              int64(zeros(0, 1)), ...
              zeros(0, 1), ...
              isLfr, bdmRemoveArray, 10, 20, false(0, 1))

            % No UFV
            test(...
              int64([10, 11, 12]'), ...
              [0, 0, 0]', ...
              isLfr, bdmRemoveArray, 10, 20, [false, false, false]' ...
              )
            % All UFV
            test(...
              int64([10, 11, 12]'), ...
              [0, 0, 0]', ...
              isLfr, bdmRemoveArray, 10, 20, [false, false, false]' ...
              )
          end
        end
      end

      % LFR
      % Test BDM=unknown/FV
      test(...
        int64([10:20]') * 1e9, ...
        [2, 2, 0, 0, 0, NaN, 0, 0, 0, 3, 3]', ...
        true, [2, 3], 1.5, 2.5, ...
        logical([1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1])' ...
        )

      % TDS
      % Test BDM=unknown/FV
      test(...
        int64([10:20]') * 1e9, ...
        [2, 2, 0, 0, 0, NaN, 0, 0, 0, 3, 3]', ...
        false, [2, 3], 1.5, 2.5, ...
        logical([1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1])' ...
        )
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Utility function for creating a simplified AsrSamplesAVoltSrm data
    % structure with identical values on every channel.
    function AsrSamplesAVoltSrm = AsSrm(v)
      assert(iscolumn(v))
      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap(...
        "bicas.proc.L1L2.AntennaSignalId", size(v, 1), 'CONSTANT', v, ...
        bicas.proc.L1L2.AntennaSignalId.ALL_ARRAY);
    end



  end    % methods(Static, Access=private)



end
