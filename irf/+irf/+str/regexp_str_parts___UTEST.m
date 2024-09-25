%
% matlab.unittest automatic test code for
% irf.str.regexp_str_parts().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2018-01-25
%
classdef regexp_str_parts___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      KEY_ANY_VALUE_REGEXP_LIST            = {'[a-z]+', '=', '.*'};   % Not whitespace before or after =.
      KEY_QUOTED_VALUE_REGEXP_LIST         = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *'};    % Require quoted value.
      KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *', '(#.*)?'};

      ES = char(zeros(1,0));    % Emty string, size 1x0.



      testCase.test_P_OK_A_OK('',     {''},       {''}',    ES, true)

      % Exact match
      % -----------
      % Non-empty regexp
      testCase.test_P_OK_A_OK('abc', {'a', 'b', 'c'},     {'a' 'b' 'c'}',    ES, true);
      % Include empty regexp.
      testCase.test_P_OK_A_OK('abc', {'a', 'b', '', 'c'}, {'a' 'b' '' 'c'}', ES, true);

      % Matching, until running out of string
      % -------------------------------------
      testCase.test_P_OK_A_exc(    'ab', {'a', 'b', 'c'}, {'a', 'b'}', ES, false);

      % First regexp matches, but second does not
      % -----------------------------------------
      % Non-empty remaining string.
      testCase.test_P_OK_A_exc(    'ac', {'a', 'b', 'c'}, {'a'}, 'c', false);
      testCase.test_P_OK_A_exc(    'ac', {'a', 'b'     }, {'a'}, 'c', false);
      % Running out of string
      testCase.test_P_OK_A_exc(    'a',  {'a', 'b', 'c'}, {'a'}, ES,  false);



      testCase.test_P_OK_A_OK('word',  {'[a-z]+'}, {'word'}', ES, true);
      testCase.test_A_exc('key = value', KEY_QUOTED_VALUE_REGEXP_LIST);

      testCase.test_P_OK_A_OK('key=value',                     KEY_ANY_VALUE_REGEXP_LIST,            {'key',                '=',         'value'             }', ES, true);
      testCase.test_P_OK_A_OK('key = "value"',                 KEY_QUOTED_VALUE_REGEXP_LIST,         {'key',               ' = ',   '"', 'value', '"', ''    }', ES, true);
      testCase.test_P_OK_A_OK('key=   "value"   ',             KEY_QUOTED_VALUE_REGEXP_LIST,         {'key',                '=   ', '"', 'value', '"', '   ' }', ES, true);
      testCase.test_P_OK_A_OK('key_name.subset   =" .-/ "   ', KEY_QUOTED_VALUE_REGEXP_LIST,         {'key_name.subset', '   =',    '"', ' .-/ ', '"', '   ' }', ES, true);
      testCase.test_P_OK_A_OK('key = "value"',                 KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, {'key',               ' = ',   '"', 'value', '"', '', ''}', ES, true);

      % Match
      testCase.test_P_OK_A_OK('key = "value"   # Comment',     KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, {'key',               ' = ',   '"', 'value', '"', '   ', '# Comment'}', ES, true);

      % Non-match
      testCase.test_P_OK_A_exc('key = "value"   % Comment',    KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, {'key',               ' = ',   '"', 'value', '"', '   ', ''}', '% Comment', false);
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)


      % Function naming convention:
      %   P = Permit non-match
      %   A = Assert match



      % Test one call.
      function test_OK(testCase, ...
          str, regexpCa, nonMatchPolicy, ...
          expSubStrCa, expRemainingStr, expIsPerfectMatch)

        [actSubStrCa, actRemainingStr, actIsPerfectMatch] = ...
          irf.str.regexp_str_parts(str, regexpCa, nonMatchPolicy);

        testCase.assertEqual(actSubStrCa,       expSubStrCa)
        testCase.assertEqual(actRemainingStr,   expRemainingStr)
        testCase.assertEqual(actIsPerfectMatch, expIsPerfectMatch)
      end



      function test_A_exc(testCase, str, regexpCa)
        testCase.assertError(...
            @() irf.str.regexp_str_parts(str, regexpCa, 'ASSERT_MATCH'), ...
            ?MException)
      end



      % Test one call for each policy, both of which should succeed.
      function test_P_OK_A_OK(testCase, ...
          str, regexpCa, ...
          expSubStrCa, expRemainingStr, expIsPerfectMatch)

        testCase.test_OK(str, regexpCa, 'PERMIT_NON_MATCH', ...
          expSubStrCa, expRemainingStr, expIsPerfectMatch)
        testCase.test_OK(str, regexpCa, 'ASSERT_MATCH', ...
          expSubStrCa, expRemainingStr, expIsPerfectMatch)
      end



      % Test one call for each policy, where PERMIT_NON_MATCH succeeds but
      % ASSERT_MATCH fails.
      function test_P_OK_A_exc(testCase, ...
          str, regexpCa, ...
          expSubStrCa, expRemainingStr, expIsPerfectMatch)

        testCase.test_OK(   str, regexpCa, 'PERMIT_NON_MATCH', ...
          expSubStrCa, expRemainingStr, expIsPerfectMatch)
        testCase.test_A_exc(str, regexpCa)
      end



  end    % methods(Access=private)



end
