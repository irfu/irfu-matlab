%
% matlab.unittest automatic test code for
% bepic.spinfit.utils.get_incrementing_array().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils_get_incrementing_array___UTEST < matlab.unittest.TestCase



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  % Technically, additional properties of testCase objects with cell array
  % default values. Test methods with arguments with the same name will be
  % called once for every element in the cell arrays.
  properties(TestParameter)
    MC = {"int64", "double", "single"}

    % Values which can be used for incrementing xRef, which in turn should not
    % affect the function return value.
    M  = {-10, -3, -2,-1, 0, 1, 2, 3, 10}
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero(T, MC)
      T.test(0, 0, 1, 0,   0, MC)
    end



    function test_0(T, MC)
      T.test(10, 10, 1, 0,   10, MC)
    end



    function test_1(T, MC)
      T.test(-10, -10, 1, 0,   -10, MC)
    end



    function test_complex_0(T, MC, M)
      T.test(10, 20, 3, 0+M*3,   [9 12 15 18]', MC)
    end

    function test_complex_1(T, MC, M)
      T.test(10, 20, 3, 1+M*3,   [10 13 16 19]', MC)
    end

    function test_complex_2(T, MC, M)
      T.test(10, 20, 3, 2+M*3,   [8 11 14 17 20]', MC)
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Run one test for input & output of a specific MATLAB CLASS ("mc").
    function test(T, xBegin, xEnd, xPeriod, xRef, expXAr, mc)
      actXAr = bepic.spinfit.utils.get_incrementing_array(...
        cast(xBegin,  mc), ...
        cast(xEnd,    mc), ...
        cast(xPeriod, mc), ...
        cast(xRef,    mc));

      T.verifyEqual( ...
        cast(actXAr, mc), ...
        cast(expXAr, mc))
    end



  end    % methods(Access=private)



end
