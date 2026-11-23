%
% matlab.unittest automatic test code for
% bepic.spinfit.ftwn.get_FTWN_fitted_center_timestamps().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ftwn_get_FTWN_fitted_center_timestamps___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_empty(T)
      FtwnTable = table();
      FtwnTable.iBegin = double.empty(0, 1);
      FtwnTable.iEnd   = double.empty(0, 1);

      actTt2000Ar = bepic.spinfit.ftwn.get_FTWN_fitted_center_timestamps( ...
        FtwnTable      = FtwnTable, ...
        tt2000Ar      = int64.empty(0, 1), ...
        spinPhaseRadAr = double.empty(0, 1));

      T.verifyEqual(actTt2000Ar, int64.empty(0, 1))
    end



    function test_0xFTWN_nonzero_samples(T)
      FtwnTable = table();
      FtwnTable.iBegin = double.empty(0, 1);
      FtwnTable.iEnd   = double.empty(0, 1);

      actTt2000Ar = bepic.spinfit.ftwn.get_FTWN_fitted_center_timestamps( ...
        FtwnTable      = FtwnTable, ...
        tt2000Ar       = int64([2; 5]), ...
        spinPhaseRadAr = [0.2; 0.5]);

      T.verifyEqual(actTt2000Ar, int64.empty(0, 1))
    end



    function test_1xFTWN(T)
      FtwnTable = table();
      FtwnTable.iBegin = [1];
      FtwnTable.iEnd   = [2];

      actTt2000Ar = bepic.spinfit.ftwn.get_FTWN_fitted_center_timestamps( ...
        FtwnTable      = FtwnTable, ...
        tt2000Ar       = int64([  2;   4]), ...
        spinPhaseRadAr =       [0.2; 0.4] * 2*pi);

      T.verifyEqual(actTt2000Ar, int64(5))
    end



    function test_2xFTWN(T)
      FtwnTable = table();
      FtwnTable.iBegin = [1; 4];
      FtwnTable.iEnd   = [2; 6];

      % NOTE: The center of the first FTWN spills over outside of the intended
      % FTWN, but the function should not check for that.
      actTt2000Ar = bepic.spinfit.ftwn.get_FTWN_fitted_center_timestamps( ...
        FtwnTable      = FtwnTable, ...
        tt2000Ar       = int64([  20;  40;  55;   60;  70;  80]), ...
        spinPhaseRadAr = [      0.2; 0.4; 1.9;  0.4; 0.6; 0.8] * 2*pi);

      T.verifyEqual(actTt2000Ar, int64([50; 65]))
    end



  end    % methods(Test)



end
