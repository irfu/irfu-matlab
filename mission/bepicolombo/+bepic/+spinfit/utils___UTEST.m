%
% matlab.unittest automatic test code for bepic.spinfit.utils (except for those
% bepic.spinfit.utils functions which have dedicated test files).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_E_ZV_epoch_to_linear(T)
      ZV_T_OFFSET = single([0, 0.25, 0.50, 0.75]);

      actTt2000Ar = bepic.spinfit.utils.E_ZV_epoch_to_linear(...
        int64.empty(0, 1), ZV_T_OFFSET);
      T.assertEqual(actTt2000Ar, int64.empty(0, 1))

      actTt2000Ar = bepic.spinfit.utils.E_ZV_epoch_to_linear(...
        int64([0]), ZV_T_OFFSET);
      T.assertEqual(actTt2000Ar, int64([0e9; 0.25e9; 0.50e9; 0.75e9]))

      actTt2000Ar = bepic.spinfit.utils.E_ZV_epoch_to_linear(...
        int64([0; 2; 5] * 1e9), ZV_T_OFFSET);
      T.assertEqual( ...
        actTt2000Ar, int64( ...
        [ ...
        0e9; 0.25e9; 0.50e9; 0.75e9; ...
        2e9; 2.25e9; 2.50e9; 2.75e9; ...
        5e9; 5.25e9; 5.50e9; 5.75e9; ...
        ]))
    end



    function test_E_ZV_spinphase_to_linear(T)
      actSpinphaseRad = bepic.spinfit.utils.E_ZV_spinphase_to_linear(...
        single.empty(0, 4));
      T.assertEqual(actSpinphaseRad, double.empty(0, 1))

      actSpinphaseRad = bepic.spinfit.utils.E_ZV_spinphase_to_linear(...
        single([0, 45, 90, 135; 180, 225, 270, 315]));
      T.assertEqual(actSpinphaseRad, double( ...
        [0; 0.25; 0.5; 0.75; 1.00; 1.25; 1.50; 1.75] * pi), "RelTol", 1e-7)
    end



    function test_E_ZV_Ev_to_linear(T)
      actSamplesAr = bepic.spinfit.utils.E_ZV_Ev_to_linear(...
        single.empty(0, 4));
      T.assertEqual(actSamplesAr, double.empty(0, 1))

      actSamplesAr = bepic.spinfit.utils.E_ZV_Ev_to_linear(...
        single([1, 2, 3, 4; 6, 7, 8, 9]));
      T.assertEqual(actSamplesAr, double([1; 2; 3; 4; 6; 7; 8; 9]))
    end



  end    % methods(Test)



end
