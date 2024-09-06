%
% matlab.unittest automatic test code for bicas.utils.SameRowsMap.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SameRowsMap___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_basic(testCase)
      % Test sequences of operations on a single SRM.
      % Exclude set_rows().

      % ==============================================
      % Adds key-values, zero rows, char keys, 'EMPTY'
      % ==============================================
      V1 = zeros(0, 1);
      V2 = ones( 0, 1, 2);

      Srm = bicas.utils.SameRowsMap('string', 0, 'EMPTY');

      testCase.assertEqual(Srm.keys(), cell(0, 1))
      testCase.assertEqual(Srm.numEntries, 0)
      testCase.assertEqual(Srm.nRows,  0)
      testCase.test_keys_values(testCase, Srm, {}, {})

      Srm.add("K1", V1)

      testCase.assertEqual(Srm.keys, {"K1"})
      testCase.assertEqual(Srm.numEntries, 1)
      testCase.test_keys_values(testCase, Srm, {"K1"}, {V1})
      testCase.assertFalse(Srm.isKey("K2"))

      Srm.add("K2", V2)

      testCase.assertEqual(Srm.numEntries, 2)
      testCase.assertTrue( Srm.isKey("K2"))

      % Test different orders. ==> Effectively testing the helper
      % function "test_keys_values()".
      testCase.test_keys_values(testCase, Srm, {"K2", "K1"}, {V2, V1})
      testCase.test_keys_values(testCase, Srm, {"K1", "K2"}, {V1, V2})

      testCase.assertEqual(Srm("K1"), V1);
      testCase.assertEqual(Srm("K2"), V2);
      testCase.assertEqual(Srm.nRows,  0)

      % =================================================
      % Zero number of constant values (test constructor)
      % =================================================
      Srm = bicas.utils.SameRowsMap('double', 3, 'CONSTANT', [1;2;3], {});
      testCase.assertEqual(Srm.nRows, 3)

      % ======================================
      % double keys, non-zero rows, 'CONSTANT'
      % ======================================
      V = [1;2;3];
      Srm = bicas.utils.SameRowsMap('double', 3, 'CONSTANT', V, {9});
      testCase.assertEqual(Srm.nRows, 3)
      testCase.assertEqual(Srm(9), V)
      testCase.test_keys_values(testCase, Srm, {9}, {V})

      % =============================================
      % Initial value has inconsistent number of rows
      % =============================================
      testCase.assertError(...
        @() (bicas.utils.SameRowsMap('double', 3, 'CONSTANT', [1;2])), ...
        ?MException)
    end



    function test_object_key_types(testCase)
      SrmDouble = bicas.utils.SameRowsMap('double', 1, 'EMPTY');
      SrmDouble.add(3.14, 123);

      SrmUint8 = bicas.utils.SameRowsMap('uint8', 1, 'EMPTY');
      SrmUint8.add(uint8(3), 123);
    end



    function test_object_key_types_error(testCase)
      SrmDouble = bicas.utils.SameRowsMap('double', 1, 'EMPTY');

      testCase.assertError(...
          @() SrmDouble.add(uint8(3), 123), ...
          ?MException)
    end



    function test_object_keys(testCase)
      Asid1 = bicas.proc.L1L2.AntennaSignalId.C.DC_V1;
      Asid2 = bicas.proc.L1L2.AntennaSignalId.C.DC_V12;
      Asid3 = bicas.proc.L1L2.AntennaSignalId.C.DC_V3;
      Srm = bicas.utils.SameRowsMap('bicas.proc.L1L2.AntennaSignalId', 1, 'EMPTY');

      Srm.add(Asid1, 1)
      Srm.add(Asid2, 1)
      testCase.assertTrue(Srm.isKey(Asid1))
      testCase.assertFalse(Srm.isKey(Asid3))
    end



    function test_set_rows(testCase)

      % Insert zero rows into zero rows.
      Srm1 = bicas.utils.SameRowsMap('string', 0, 'EMPTY');
      Srm1.add("K2", zeros(0,0))
      Srm2 = bicas.utils.SameRowsMap('string', 0, 'EMPTY');
      Srm2.add("K2", zeros(0,0))
      Srm1.set_rows(Srm2, zeros(0,1))
      testCase.assertEqual(Srm1("K2"), zeros(0,0))

      % Insert zero rows into non-zero rows.
      Srm1 = bicas.utils.SameRowsMap('string', 3, 'EMPTY');
      Srm1.add("K2", zeros(3,0))
      Srm2 = bicas.utils.SameRowsMap('string', 0, 'EMPTY');
      Srm2.add("K2", zeros(0,0))
      Srm1.set_rows(Srm2, zeros(0,1))
      testCase.assertEqual(Srm1("K2"), zeros(3,0))

      % Preserve type
      Srm1 = bicas.utils.SameRowsMap('string', 4, 'EMPTY');
      Srm1.add("K2", int16([1;2;3;4]))
      Srm2 = bicas.utils.SameRowsMap('string', 2, 'EMPTY');
      Srm2.add("K2", int16([[8;9]]))
      Srm1.set_rows(Srm2, [2;3])
      testCase.assertEqual(Srm1("K2"), int16([1;8;9;4]))



      % =======================
      % Higher-dimensional case
      % =======================
      % IMPLEMENTATION NOTE: Simpler to just define one example variable
      % (of higher dimensionality) and then extract smaller variables from
      % it.
      V1(:, :, 1) = [...
        111, 121, 131; ...
        211, 221, 231; ...
        311, 321, 331; ...
        411, 421, 431 ...
        ];
      V1(:, :, 2) = [...
        112, 122, 132; ...
        212, 222, 232; ...
        312, 322, 332; ...
        412, 422, 432 ...
        ];

      Srm1 = bicas.utils.SameRowsMap('string', 4, 'CONSTANT', V1,            {"K"});
      Srm2 = bicas.utils.SameRowsMap('string', 2, 'CONSTANT', V1(1:2, :, :), {"K"});
      Srm1.set_rows(Srm2, [3;2])    % NOTE: Decrementing indices.

      V2              = V1;
      V2([3,2], :, :) = V1(1:2, :, :);
      testCase.assertEqual(Srm1("K"), V2)



      % =======
      % Illegal
      % =======
      % Incompatible types
      Srm1 = bicas.utils.SameRowsMap('string', 2, 'CONSTANT', int8([1;2]), {"K"});
      Srm2 = bicas.utils.SameRowsMap('string', 1, 'CONSTANT', [9],         {"K"});
      testCase.assertError(...
        @() (Srm1.set_rows(Srm2, [1])), ...
        ?MException)
      % Incompatible array sizes
      Srm1 = bicas.utils.SameRowsMap('string', 2, 'CONSTANT', [1 2;3 4], {"K"});
      Srm2 = bicas.utils.SameRowsMap('string', 1, 'CONSTANT', [9],       {"K"});
      testCase.assertError(...
        @() (Srm1.set_rows(Srm2, [1])), ...
        ?MException)
      % Different sets of keys.
      Srm1 = bicas.utils.SameRowsMap('string', 2, 'CONSTANT', [1;2], {"K1", "K2a"});
      Srm2 = bicas.utils.SameRowsMap('string', 1, 'CONSTANT', [9],   {"K1", "K2b"});
      testCase.assertError(...
        @() (Srm1.set_rows(Srm2, [1])), ...
        ?MException)
    end



    % NOTE: Does not test method per se.
    function test_equal(testCase)

      % NOTE: Having no keys with different KeyType means that the key
      % values do not distinguish the KeyTypes.
      Srm1a = bicas.utils.SameRowsMap('string', 1, 'EMPTY');
      Srm1b = bicas.utils.SameRowsMap('string', 1, 'EMPTY');
      Srm2a = bicas.utils.SameRowsMap('double', 1, 'EMPTY');
      Srm2b = bicas.utils.SameRowsMap('double', 1, 'EMPTY');
      Srm3  = bicas.utils.SameRowsMap('string', 2, 'EMPTY');

      testCase.assertTrue( Srm1a == Srm1b)
      testCase.assertFalse(Srm1a == Srm2a)
      testCase.assertTrue( Srm2a == Srm2b)
      testCase.assertFalse(Srm1a == Srm3 )

      Srm1a = bicas.utils.SameRowsMap('string', 1, 'CONSTANT', [1],   {"K1", "K2"});
      Srm1b = bicas.utils.SameRowsMap('string', 1, 'CONSTANT', [1],   {"K1", "K2"});
      Srm2a = bicas.utils.SameRowsMap('double', 1, 'CONSTANT', [9],   {1, 2});
      Srm2b = bicas.utils.SameRowsMap('double', 1, 'CONSTANT', [9],   {1, 2});
      Srm3  = bicas.utils.SameRowsMap('string', 1, 'CONSTANT', [9],   {"K1", "K2"});
      Srm4  = bicas.utils.SameRowsMap('string', 2, 'CONSTANT', [1;2], {"K1", "K2"});

      testCase.assertTrue( Srm1a == Srm1b)
      testCase.assertTrue( Srm2a == Srm2b)
      testCase.assertFalse(Srm1a == Srm3 )
      testCase.assertFalse(Srm1a == Srm4)
      testCase.assertFalse(Srm2a == Srm3)
      testCase.assertFalse(Srm2a == Srm4)

      % NaN, different key order.
      Srm1a = bicas.utils.SameRowsMap('string', 1, 'EMPTY');
      Srm1a.add("K1", [NaN])
      Srm1a.add("K2", [2])
      Srm1b = bicas.utils.SameRowsMap('string', 1, 'EMPTY');
      Srm1b.add("K2", [2])
      Srm1b.add("K1", [NaN])

      testCase.assertTrue( Srm1a == Srm1b)
      testCase.assertTrue( Srm1a == Srm1b)
      testCase.assertFalse(Srm1a == Srm3 )

      % Different MATLAB classes.
      Srm1  = bicas.utils.SameRowsMap('string', 1, 'EMPTY');
      Srm1.add("K1", [1])
      Srm3  = bicas.utils.SameRowsMap('string', 1, 'EMPTY');
      Srm3.add("K1", int8([1]))
      testCase.assertFalse(Srm1a == Srm3 )
    end



    function test_display(testCase)
      Srm = bicas.utils.SameRowsMap('string', 0, 'CONSTANT', uint8(ones(0,2)), {"K1"});
      disp(Srm)

      Srm = bicas.utils.SameRowsMap('double', 1, 'CONSTANT', [9], {1});
      disp(Srm)

      Srm = bicas.utils.SameRowsMap('string', 3, 'CONSTANT', [1,2;3,4;5,6], {"K1", "K2"});
      disp(Srm)

      Asid1 = bicas.proc.L1L2.AntennaSignalId.C.DC_V1;
      Asid2 = bicas.proc.L1L2.AntennaSignalId.C.DC_V12;
      Srm = bicas.utils.SameRowsMap('string', 1, 'CONSTANT', [Asid1, Asid2], {"K1", "K2"});
      disp(Srm)
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Helper function. Test that SRM returns correct values for methods
    % .keys() and .values(). Takes into account that keys/values may be in
    % a different order than specified.
    function test_keys_values(testCase, Srm, expKeysCa, expValuesCa)
      % NOTE: Implementation should permit both numbers and strings as
      % keys, mixed.

      assert(isa(Srm, 'bicas.utils.SameRowsMap'))

      actKeysCa   = Srm.keys();
      actValuesCa = Srm.values();

      testCase.assertTrue(bicas.utils.object_sets_isequaln(actKeysCa, expKeysCa))

      for iAct = 1:numel(actKeysCa)
        for iExp = 1:numel(expKeysCa)
          if isequal(actKeysCa{iAct}, expKeysCa{iExp})
            actValue = actValuesCa{iAct};
            expValue = expValuesCa{iExp};
            testCase.assertEqual(actValue, expValue)
          end
        end
      end

    end



  end    % methods(Static, Access=private)



end
