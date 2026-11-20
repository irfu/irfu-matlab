%
% matlab.unittest automatic test code for
% bepic.spinfit.fit_Efield_SAFW().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_fit_Efield_SAFW___UTEST < matlab.unittest.TestCase
  % PROPOSAL: More tests.
  % PROPOSAL: Test data gaps creating too small fit windows.



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    N_COEFF_OFFICIAL = bepic.spinfit.N_COEFF_OFFICIAL_EFIELD_FIT
    N_COEFF_IRFU     = bepic.spinfit.N_COEFF_IRFU_EFIELD_FIT
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_samples(T)
      ActT = bepic.spinfit.fit_Efield_SAFW( ...
        tt2000Ar                = int64.empty(0, 1), ...
        spinPhaseRadAr          = double.empty(0, 1), ...
        samplesAr               = double.empty(0, 1), ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF_OFFICIAL);

      ExpT = table();
      ExpT.A                     = double.empty(0, T.N_COEFF_OFFICIAL);
      ExpT.AIrfu                 = double.empty(0, T.N_COEFF_IRFU);
      ExpT.fitWindowCenterTt2000 = int64.empty( 0, 1);

      T.verifyEqual(ActT, ExpT)
    end



    function test_one_sample(T)
      % No output
      ActT = bepic.spinfit.fit_Efield_SAFW( ...
        tt2000Ar                = int64([10]'), ...
        spinPhaseRadAr          = [0.5]' * 2*pi, ...
        samplesAr               = [3], ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF_OFFICIAL);

      ExpT = table();
      ExpT.A                     = double.empty(0, T.N_COEFF_OFFICIAL);
      ExpT.AIrfu                 = double.empty(0, T.N_COEFF_IRFU);
      ExpT.fitWindowCenterTt2000 = int64.empty( 0, 1);

      T.verifyEqual(ActT, ExpT)
    end



    % Generates NaN return value.
    function test_two_samples_per_FTWN(T)
      actT = bepic.spinfit.fit_Efield_SAFW( ...
        tt2000Ar                = int64(    [ 01   06   11   16]'), ...
        spinPhaseRadAr          = wrapTo2Pi([0.1, 0.6, 0.1, 0.6]' * 2*pi), ...
        samplesAr               =           [  1; 2; 3; 4], ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF_OFFICIAL);

      T.verifyEqual(actT.fitWindowCenterTt2000, int64([5; 15]))
      T.verifyEqual(actT.A,                     NaN(2, T.N_COEFF_OFFICIAL))
      T.verifyEqual(actT.AIrfu,                 NaN(2, T.N_COEFF_IRFU))
    end



    % Generates finite+NaN return value.
    function test_finite_RV(T)
      actT = bepic.spinfit.fit_Efield_SAFW( ...
        tt2000Ar                = int64([ 0:21 ]'), ...
        spinPhaseRadAr          = mod(  [0.0 : 0.1 : 2.1]' * 2*pi, 2*pi), ...
        samplesAr               =       [  0:21 ]', ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = max(T.N_COEFF_OFFICIAL, T.N_COEFF_IRFU));

      T.verifyEqual(actT.fitWindowCenterTt2000, int64([5; 15; 25]))
      T.verifyEqual(isfinite(actT.A),           [true(2, T.N_COEFF_OFFICIAL); false(1, T.N_COEFF_OFFICIAL)])
      T.verifyEqual(isfinite(actT.AIrfu),       [true(2, T.N_COEFF_IRFU);     false(1, T.N_COEFF_IRFU    )])
    end



    % Generates finite+NaN return value + data gap.
    function test_data_gap_finite_RV(T)
      actT = bepic.spinfit.fit_Efield_SAFW( ...
        tt2000Ar                = int64([ [0:9] 1000+[10:21] ]'), ...
        spinPhaseRadAr          = mod(  [0.0 : 0.1 : 2.1]' * 2*pi, 2*pi), ...
        samplesAr               =       [  0:21 ]', ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = max(T.N_COEFF_OFFICIAL, T.N_COEFF_IRFU));

      T.verifyEqual(actT.fitWindowCenterTt2000, int64([5; 1015; 1025]))
      T.verifyEqual(isfinite(actT.A),           [true(2, T.N_COEFF_OFFICIAL); false(1, T.N_COEFF_OFFICIAL)])
      T.verifyEqual(isfinite(actT.AIrfu),       [true(2, T.N_COEFF_IRFU);     false(1, T.N_COEFF_IRFU    )])
    end



  end    % methods(Test)



end
