%
% matlab.unittest automatic test code for bicas.utils.ZvMap.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ZvMap___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Test adding custom-made class as value.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_add_set_get_1(T)
      Zvm = bicas.utils.ZvMap(0);

      % Illegal number of value rows.
      T.assertError(...
        @() Zvm.add(double(3), ones(2,2)), ...
        ?MException)

      T.test_add_set_get_OK(Zvm, "KEY",     uint8(ones(0, 2, 3)), int16(ones(0, 3, 2)))

      % Distinguish between numerically equal keys.
      T.test_add_set_get_OK(Zvm, int8(3),   ones(0, 0, 1), ones(0, 1, 0))
      T.test_add_set_get_OK(Zvm, uint8(3),  ones(0, 2, 3), ones(0, 3, 2))
      T.test_add_set_get_OK(Zvm, uint16(3), ones(0, 4, 5), ones(0, 5, 4))
      T.test_add_set_get_OK(Zvm, double(3), ones(0, 6, 7), ones(0, 7, 6))

      % Pre-existing key.
      T.assertError(...
        @() Zvm.add(uint16(3), ones(0, 2)), ...
        ?MException)
    end



    function test_add_set_get_2(T)
      Zvm = bicas.utils.ZvMap(3);

      T.test_add_set_get_OK(Zvm, uint64(3), ...
        [1;2;3], [3;2;1])
      T.test_add_set_get_OK(Zvm, uint16(3), ...
        [4 5; 6 7; 8 9], [8 9; 6 7; 4 5])
    end



    function test_nEntries_keyCa(T)
      Zvm = bicas.utils.ZvMap(4);

      T.assertEqual(Zvm.nEntries, 0)
      T.assertEqual(Zvm.keyCa,    cell(0, 1))

      Zvm.add(3, ones(4,2))
      T.assertEqual(Zvm.nEntries, 1)
      T.assertEqual(Zvm.keyCa,    {3})

      Zvm.add("ABC", [2;3;4;5])
      T.assertEqual(Zvm.nEntries, 2)
      T.assertEqual(Zvm.keyCa,    {3; "ABC"})
    end



    function test_equality_empty(T)
      Zvm1a = bicas.utils.ZvMap(3);
      Zvm1b = bicas.utils.ZvMap(3);
      Zvm2  = bicas.utils.ZvMap(4);

      T.test_equality_helper(Zvm1a, Zvm1b, Zvm2)
    end



    function test_equality_nonempty_1(T)
      Zvm1a = bicas.utils.ZvMap(1);
      Zvm1b = bicas.utils.ZvMap(1);
      Zvm2  = bicas.utils.ZvMap(1);
      Zvm1a.add("KEY", [1, 2, 3])
      Zvm1b.add("KEY", [1, 2, 3])
      Zvm2.add( "KEY", [1, 2, 4])

      T.test_equality_helper(Zvm1a, Zvm1b, Zvm2)
    end



    function test_equality_nonempty_2(T)
      Zvm1a = bicas.utils.ZvMap(1);
      Zvm1b = bicas.utils.ZvMap(1);
      Zvm2  = bicas.utils.ZvMap(1);

      VALUE = [1, 2, 3];
      Zvm1a.add("KEY",  VALUE)
      Zvm1b.add("KEY",  VALUE)
      Zvm2.add( "KEY2", VALUE)

      T.test_equality_helper(Zvm1a, Zvm1b, Zvm2)
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % NOTE: Also indirectly tests whether the object is actually modified
    % without the variable needing to be assigned, i.e. whether it is a
    % (working) handle class.
    function test_add_set_get_OK(T, Zvm, key, value1, value2)
      assert(~isequaln(value1, value2))

      % add(), get()
      Zvm.add(key, value1)
      actValue1 = Zvm.get(key);
      T.assertEqual(actValue1, value1)

      % set(), get()
      Zvm.set(key, value2)
      actValue2 = Zvm.get(key);
      T.assertEqual(actValue2, value2)
    end



    % Test both equality and inequality.
    function test_equality_helper(T, Zvm1a, Zvm1b, Zvm2)
      T.assertTrue(isequaln( Zvm1a, Zvm1b))
      T.assertEqual(         Zvm1a, Zvm1b)

      T.assertFalse(isequaln(Zvm1a, Zvm2))
      T.assertNotEqual(      Zvm1a, Zvm2)
    end



  end    % methods(Access=private)



end
