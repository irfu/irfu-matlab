%
% matlab.unittest automatic test code for bicas.proc.L2L3.ext.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ext___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_calc_EFIELD_SCPOT(T)

      % % Arbitrary number output variables.
      % function test(inputsCa, expOutputsCa)
      %   % Pre-allocate correct size for later assignment via function.
      %   actOutputs = cell(size(expOutputsCa));
      %
      %   [actOutputs{:}] = FUNCTION_TO_TEST(inputsCa{:});
      %   T.assertEqual(actOutputs, expOutputsCa)
      % end
      %
      % function test_exc(varargin)
      %   T.assertError(...
      %       @() FUNCTION_TO_TEST(varargin{:}), ...
      %       ?MException)
      % end
      %===================================================================

      % Function input
      tt2000  = int64([10:14]');
      VDC_Fpa = bicas.utils.FPArray(single(zeros(numel(tt2000), 3)));
      EDC_Fpa = bicas.utils.FPArray(single(zeros(numel(tt2000), 3)));

      % =========================
      % EXCD/solo.vdccal() output
      % =========================
      VdccalRv.DCE_SRF_out = irf.ts_vec_xyz( ...
        EpochTT(tt2000), ...
        [ ...
        NaN,   1,   2; ...
        NaN,   1,   2; ...
        NaN, NaN,   2; ...
        NaN,   1, NaN; ...
        NaN, NaN, NaN; ...
        ] ...
        );
      VdccalRv.DCE_SRF_out.units            = 'mV/m';
      VdccalRv.DCE_SRF_out.coordinateSystem = 'SRF';

      VdccalRv.PSP_out = TSeries( ...
        EpochTT(tt2000), ...
        [ ...
        1; ...
        2; ...
        NaN;
        4;
        5;
        ]);
      VdccalRv.PSP_out.units = 'V';

      VdccalRv.ScPot_out = TSeries( ...
        EpochTT(tt2000), ...
        [ ...
        1; ...
        NaN; ...
        3; ...
        4; ...
        5; ...
        ]);
      VdccalRv.ScPot_out.units = 'V';

      VdccalRv.codeVerStr = '2026-05-19T14:15:00';
      VdccalRv.matVerStr  = 'd23K123_20230707.mat';
      Psp2neRv = {};
      Excd = bicas.proc.L2L3.ExternalCodeTest(VdccalRv, Psp2neRv);

      % ================
      % CALL TESTED CODE
      % ================
      actR = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT(...
        Excd, tt2000=tt2000, VDC_Fpa=VDC_Fpa, EDC_Fpa=EDC_Fpa);

      assert(isequaln(...
        actR.EdcSrfTs.data, ...
        [
        NaN,   1,   2; ...
        NaN,   1,   2; ...
        NaN, NaN, NaN; ...
        NaN, NaN, NaN; ...
        NaN, NaN, NaN; ...
        ] ...
      ))
    end



  end    % methods(Test)



end
