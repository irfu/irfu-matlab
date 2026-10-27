%
% matlab.unittest automatic test code for bepic.spinfit_fit_array.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_fit_array___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Reuse data somehow, and only modify fit windows to generate
  %           different cases?



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  % Technically, additional properties of testCase objects with cell array
  % default values. Test methods with arguments with the same name will be
  % called once for every element in the cell arrays.
  properties(TestParameter)
    % nMinSamplesPerFtwn for when the value should not matter.
    N_MIN_SAMPLES_PER_FTWN = {0, 1, 2, 10}
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_fit_windows_zero_samples(T, N_MIN_SAMPLES_PER_FTWN)
      FtwnTable = table();
      FtwnTable.beginTt2000 = int64.empty(0, 1);
      FtwnTable.endTt2000   = int64.empty(0, 1);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar           = int64.empty(0, 1), ...
        samplesAr          = double.empty(0, 1), ...
        spinPhaseRadAr     = double.empty(0, 1), ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = @(C, x) (C(1)*x + C(2)), ...
        fittingStartValue  = [0, 0], ...
        nMinSamplesPerFtwn = N_MIN_SAMPLES_PER_FTWN);

      ExpT = table();
      ExpT.coefficients = double.empty(0, 2);

      T.verifyEqual(ActT, ExpT)
    end



    function test_one_fit_window_zero_samples(T, N_MIN_SAMPLES_PER_FTWN)
      FtwnTable = table();
      FtwnTable.beginTt2000 = int64(10);
      FtwnTable.endTt2000   = int64(20);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar           = int64.empty( 0, 1), ...
        spinPhaseRadAr     = double.empty(0, 1), ...
        samplesAr          = double.empty(0, 1), ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = @(C, x) (C(1)*x + C(2)), ...
        fittingStartValue  = [0, 0], ...
        nMinSamplesPerFtwn = N_MIN_SAMPLES_PER_FTWN);

      ExpT = table();
      ExpT.coefficients = [NaN, NaN];

      T.verifyEqual(ActT, ExpT)
    end



    % More than zero samples, but they are too few.
    function test_one_fit_window_too_few_samples(T, N_MIN_SAMPLES_PER_FTWN)
      FtwnTable = table();
      FtwnTable.beginTt2000 = int64(14);
      FtwnTable.endTt2000   = int64(16);
      % ==> One sample in one fit window.

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar           = int64([5 10 15 20 25])', ...
        spinPhaseRadAr     = [0.1 0.2 0.3 0.4 0.5]', ...
        samplesAr          = [  1   2   3   4   5]', ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = @(C, x) (C(1)*x + C(2)), ...
        fittingStartValue  = [0, 0], ...
        nMinSamplesPerFtwn = N_MIN_SAMPLES_PER_FTWN);

      ExpT = table();
      ExpT.coefficients = [NaN, NaN];

      T.verifyEqual(ActT, ExpT)
    end



    % Some samples are NaN, but the remaining ones are enough.
    % Test whether nMinSamplesPerFtwn makes a difference.
    function test_one_fit_window_NaN_samples(T)

      function ActT = test(nMinSamplesPerFtwn)
        FtwnTable = table();
        FtwnTable.beginTt2000 = int64(5);
        FtwnTable.endTt2000   = int64(25);

        ActT = bepic.spinfit.fit_array( ...
          tt2000Ar           = int64([5 10 15 20 25])', ...
          spinPhaseRadAr     = [0.1 0.2 0.3 0.4 0.5]', ...
          samplesAr          = [ 31 NaN  33 NaN  55]', ...
          FtwnTable          = FtwnTable, ...
          fittingFunction    = @(C, x) (C(1)*x + C(2)), ...
          fittingStartValue  = [0, 0], ...
          nMinSamplesPerFtwn = nMinSamplesPerFtwn);
      end

      ActT = test(2);
      ExpT = table();
      ExpT.coefficients = [10, 30];

      T.verifyEqual(ActT, ExpT, "AbsTol", 1e-5)



      ActT = test(3);
      ExpT = table();
      ExpT.coefficients = [NaN, NaN];

      T.verifyEqual(ActT, ExpT, "AbsTol", 1e-5)
    end



    function test_one_fit_window_enough_samples(T)
      FtwnTable = table();
      FtwnTable.beginTt2000 = int64(10);
      FtwnTable.endTt2000   = int64(20);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar           = int64([5 10 15 20 25])', ...
        spinPhaseRadAr     = [0.1 0.2 0.3 0.4 0.5]', ...
        samplesAr          = [ 31  32  33  34  35]', ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = @(C, x) (C(1)*x + C(2)), ...
        fittingStartValue  = [0, 0], ...
        nMinSamplesPerFtwn = 2);

      ExpT = table();
      ExpT.coefficients = [10, 30];

      T.verifyEqual(ActT, ExpT, "AbsTol", 1e-5)
    end



    function test_two_fit_windows_enough_samples(T)
      FtwnTable = table();
      FtwnTable.beginTt2000 = int64([5  15]');
      FtwnTable.endTt2000   = int64([15 25]');

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar                = int64([5 10 15 20 25])', ...
        spinPhaseRadAr     = [0.1 0.2 0.3 0.4 0.5]', ...
        samplesAr          = [ 31  32  53  54  75]', ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = @(C, x) (C(1)*x + C(2)), ...
        fittingStartValue  = [0, 0], ...
        nMinSamplesPerFtwn = 2);

      ExpT = table();
      ExpT.coefficients = [10, 30; 10, 50];

      T.verifyEqual(ActT, ExpT, "AbsTol", 1e-5)
    end



  end    % methods(Test)



end
