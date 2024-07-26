%
% matlab.unittest automatic test code for bicas.utils.
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



    function test_get_paths(testCase)
      path = bicas.utils.get_BICAS_root_dir();
      irf.assert.dir_exists(path)

      path = bicas.utils.get_SWD_file();
      irf.assert.file_exists(path)

      path = bicas.utils.get_BICAS_config_dir();
      irf.assert.dir_exists(path)

      ~ = bicas.utils.get_BICAS_default_config_file();
      % NOTE: Can not test whether the default config file exists or not, since
      % it depends on the developer setup, and whether running in CI or not. The
      % version delivered to ROC has no default config file (ROC has to specify
      % their own).
    end



    function test_object_sets_isequaln(testCase)

      % One output variable.
      function test(keysCa1, keysCa2, expEqual)
        actEqual = bicas.utils.object_sets_isequaln(keysCa1, keysCa2);
        testCase.assertEqual(actEqual, expEqual)
      end

      test({}, {}, true)
      test({ 1 }, { 1 }, true)
      test({'1'}, {'1'}, true)
      test({ 1 }, { 2 }, false)
      test({'1'}, {'2'}, false)

      test({'asd', 1}, {'asd', 1}, true)
      test({'asd', 1}, {'asd', 2}, false)
      test({'asd', 1}, {'ASD', 1}, false)
    end



  end    % methods(Test)



end
