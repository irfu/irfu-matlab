%
% matlab.unittest automatic test code for bicas.swm.get_SWML().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_SWML___UTEST < matlab.unittest.TestCase



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  % Technically, additional properties of testCase objects with cell array
  % default values. Test methods with arguments with the same name will be
  % called once for every element in the cell arrays.
  properties(TestParameter)
    L1L2_ENABLED      = {false; true}
    L2L2_ENABLED      = {false; true}
    L2L3_SURV_ENABLED = {false; true}
    L2L3_SBMx_ENABLED = {false; true}
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Only test that SWML can be generated (without crashing) for multiple
    % combinations of settings.
    function test_combinations_of_settings(T, ...
        L1L2_ENABLED, L2L2_ENABLED, L2L3_SURV_ENABLED, L2L3_SBMx_ENABLED)

      function test(l1l2Enabled, l2l2Enabled, l2l3Enabled, l2l3SbmxEnabled)
        % Tests (1) non-crash, (2) class of return value.
        Bso = T.get_BSO(...
          l1l2Enabled, l2l2Enabled, l2l3Enabled, l2l3SbmxEnabled);

        Swml = bicas.swm.get_SWML(Bso);

        T.verifyClass(Swml, 'bicas.swm.SoftwareModeList')
      end

      %===================================================================

      test(L1L2_ENABLED, L2L2_ENABLED, L2L3_SURV_ENABLED, L2L3_SBMx_ENABLED)
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function Bso = get_BSO(l1l2Enabled, l2l2Enabled, l2l3SurvEnabled, l2l3SbmxEnabled)
      Bso = bicas.create_default_BSO();

      Bso.override_value('SWM.L1-L2_ENABLED',         l1l2Enabled,     'test')
      Bso.override_value('SWM.L2-L2_CWF-DSR_ENABLED', l2l2Enabled,     'test')
      Bso.override_value('SWM.L2-L3_SURV_ENABLED',    l2l3SurvEnabled, 'test')
      Bso.override_value('SWM.L2-L3_SBMx_ENABLED',    l2l3SbmxEnabled, 'test')

      Bso.make_read_only()
    end



  end    % methods(Static, Access=private)



end
