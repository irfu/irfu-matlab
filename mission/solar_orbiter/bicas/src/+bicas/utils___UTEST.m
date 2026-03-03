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

      [~] = bicas.utils.get_BICAS_default_config_file();
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



    function test_dict_lookup___no_error(testCase)
      % TODO: NaN key

      % Test values meant to be type cast later.
      DICT_KEY_AR   = [10, 11];
      DICT_VALUE_AR = [20, 21];
      KEY_AR        = [10 11; 10 11; 11 10];
      VALUE_AR      = [20 21; 20 21; 21 20];

      % BICAS's primary use case.
      test(...
        uint8(DICT_KEY_AR), ...
        uint8(DICT_VALUE_AR), ...
        uint8(KEY_AR), ...
        uint8(VALUE_AR));

      % Different key/value types.
      test(...
        uint8(DICT_KEY_AR), ...
        int16(DICT_VALUE_AR), ...
        uint8(KEY_AR), ...
        int16(VALUE_AR));

      test(...
        uint8(DICT_KEY_AR), ...
        uint16(DICT_VALUE_AR), ...
        uint8(KEY_AR), ...
        uint16(VALUE_AR));

      % Test non-number MATLAB classs.
      test(...
        uint8([10 11]), ...
        string(["abc" "defgh"]), ...
        uint8([10 11; 11 10; 10 10]), ...
        ["abc" "defgh"; "defgh" "abc"; "abc" "abc"]);
      test(...
        string(["abc" "defgh"]), ...
        uint8([10 11]), ...
        ["abc" "defgh"; "defgh" "abc"; "abc" "abc"], ...
        uint8([10 11; 11 10; 10 10]));

      % *NOT* 1-to-1 dictionary (not bijective)
      test(...
        uint8([10 11 12]), ...
        int16([20 21 21]), ...
        uint8([10 11 12; 12 11 10]), ...
        int16([20 21 21; 21 21 20]));

      % Empty keys (zero elements).
      for sizeCa = {[0,0], [3,0], [0,3]}
        % Empty key/value array. Non-empty dictionary.
        test(...
          uint8(DICT_KEY_AR), ...
          uint8(DICT_VALUE_AR), ...
          uint8(zeros(sizeCa{1})), ...
          uint8(zeros(sizeCa{1})));

        % Define key/value MATLAB classes. Empty dictionary.
        test(...
          uint8([]), ...
          int16([]), ...
          uint8(zeros(sizeCa{1})), ...
          int16(zeros(sizeCa{1})));
      end


      % Test successful call.
      function test(dictKeyAr, dictValueAr, keyAr, expValueAr)
        dict = dictionary(dictKeyAr, dictValueAr);

        actValueAr = bicas.utils.dict_lookup(dict, keyAr);

        testCase.assertEqual(actValueAr, expValueAr)
      end
    end



    function test_dict_lookup___error(testCase)
      % Undefined key/value MATLAB classes. Zero entries.
      % Non-zero keys.
      dict = dictionary;
      test_exc(dict, uint8([10 11]))

      % Undefined key/value MATLAB classes. Zero entries.
      % Zero keys.
      dict = dictionary;
      test_exc(dict, uint8([]))

      % Key not in non-empty dictionary.
      dict = dictionary(uint8([10 11]), int16([20 21]));
      test_exc(dict, uint8([99]))

      function test_exc(dict, keyAr)
        testCase.assertError(...
          @() bicas.utils.dict_lookup(dict, keyAr), ...
          ?MException)
      end
    end



  end    % methods(Test)



end
