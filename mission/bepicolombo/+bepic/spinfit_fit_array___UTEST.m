%
% matlab.unittest automatic test code for bepic.spinfit_fit_array.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_fit_array___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Reuse data somehow, and only modify fit windows to generate
  %           different cases?



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_fit_windows_zero_samples(T)
      FitWindowsTable = table();
      FitWindowsTable.beginTt2000 = int64.empty(0, 1);
      FitWindowsTable.endTt2000   = int64.empty(0, 1);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar         =int64.empty(0, 1), ...
        samplesAr        =double.empty(0, 1), ...
        spinPhaseRadAr   =double.empty(0, 1), ...
        FitWindowsTable  =FitWindowsTable, ...
        fittingFunction  =@(C, x) (C(1)*x + C(2)), ...
        fittingStartValue=[0, 0]);

      ExpT = table();
      ExpT.coefficients = double.empty(0, 2);

      T.verifyEqual(ActT, ExpT)
    end



    function test_one_fit_window_zero_samples(T)
      FitWindowsTable = table();
      FitWindowsTable.beginTt2000 = int64(10);
      FitWindowsTable.endTt2000   = int64(20);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar         =int64.empty( 0, 1), ...
        spinPhaseRadAr   =double.empty(0, 1), ...
        samplesAr        =double.empty(0, 1), ...
        FitWindowsTable  =FitWindowsTable, ...
        fittingFunction  =@(C, x) (C(1)*x + C(2)), ...
        fittingStartValue=[0, 0]);

      ExpT = table();
      ExpT.coefficients = [NaN, NaN];

      T.verifyEqual(ActT, ExpT)
    end



    % More that zero samples, but they are too few.
    function test_one_fit_window_too_few_samples(T)
      FitWindowsTable = table();
      FitWindowsTable.beginTt2000 = int64(14);
      FitWindowsTable.endTt2000   = int64(16);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar         =int64([5 10 15 20 25])', ...
        spinPhaseRadAr   =[0.1 0.2 0.3 0.4 0.5]', ...
        samplesAr        =[  1   2   3   4   5]', ...
        FitWindowsTable  =FitWindowsTable, ...
        fittingFunction  =@(C, x) (C(1)*x + C(2)), ...
        fittingStartValue=[0, 0]);

      ExpT = table();
      ExpT.coefficients = [NaN, NaN];

      T.verifyEqual(ActT, ExpT)
    end



    % Some samples are NaN, but the remaining ones are enough.
    function test_one_fit_window_NaN_samples(T)
      FitWindowsTable = table();
      FitWindowsTable.beginTt2000 = int64(5);
      FitWindowsTable.endTt2000   = int64(25);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar         =int64([5 10 15 20 25])', ...
        spinPhaseRadAr   =[0.1 0.2 0.3 0.4 0.5]', ...
        samplesAr        =[ 31 NaN  33 NaN  55]', ...
        FitWindowsTable  =FitWindowsTable, ...
        fittingFunction  =@(C, x) (C(1)*x + C(2)), ...
        fittingStartValue=[0, 0]);

      ExpT = table();
      ExpT.coefficients = [10, 30];

      T.verifyEqual(ActT, ExpT, "AbsTol", 1e-5)
    end



    function test_one_fit_window_enough_samples(T)
      FitWindowsTable = table();
      FitWindowsTable.beginTt2000 = int64(10);
      FitWindowsTable.endTt2000   = int64(20);

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar         =int64([5 10 15 20 25])', ...
        spinPhaseRadAr   =[0.1 0.2 0.3 0.4 0.5]', ...
        samplesAr        =[ 31  32  33  34  35]', ...
        FitWindowsTable  =FitWindowsTable, ...
        fittingFunction  =@(C, x) (C(1)*x + C(2)), ...
        fittingStartValue=[0, 0]);

      ExpT = table();
      ExpT.coefficients = [10, 30];

      T.verifyEqual(ActT, ExpT, "AbsTol", 1e-5)
    end



    function test_two_fit_windows_enough_samples(T)
      FitWindowsTable = table();
      FitWindowsTable.beginTt2000 = int64([5  15]');
      FitWindowsTable.endTt2000   = int64([15 25]');

      ActT = bepic.spinfit.fit_array( ...
        tt2000Ar         =int64([5 10 15 20 25])', ...
        spinPhaseRadAr   =[0.1 0.2 0.3 0.4 0.5]', ...
        samplesAr        =[ 31  32  53  54  75]', ...
        FitWindowsTable  =FitWindowsTable, ...
        fittingFunction  =@(C, x) (C(1)*x + C(2)), ...
        fittingStartValue=[0, 0]);

      ExpT = table();
      ExpT.coefficients = [10, 30; 10, 50];

      T.verifyEqual(ActT, ExpT, "AbsTol", 1e-5)
    end



  end    % methods(Test)



end
