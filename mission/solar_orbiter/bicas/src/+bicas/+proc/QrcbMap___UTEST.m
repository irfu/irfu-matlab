%
% matlab.unittest automatic test code for bicas.proc.QrcbMap.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcbMap___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_qrcidAr(testCase)
      QrcbMap = bicas.proc.QrcbMap(2);
      testCase.assertEqual(QrcbMap.qrcidAr, string.empty(0, 1))

      % IMPLEMENTATION NOTE: Deliberately add QRCIDs in non-alphanumeric order.
      QrcbMap.add("QRCID2", logical([0; 1]))
      testCase.assertEqual(QrcbMap.qrcidAr, ["QRCID2"])

      % Deliberately assert alphanumeric order for return value.
      QrcbMap.add("QRCID1", logical([1; 0]))
      testCase.assertEqual(QrcbMap.qrcidAr, ["QRCID1"; "QRCID2"])
    end



    function test_add_set_get_has_QRCID(testCase)
      QRCB_AR = logical([0; 1]);
      QrcbMap = bicas.proc.QrcbMap(2);

      testCase.assertError(...
          @() QrcbMap.set("QRCID1", ~QRCB_AR), ...
          ?MException)

      testCase.assertFalse(QrcbMap.has_QRCID("QRCID1"))

      QrcbMap.add("QRCID1", QRCB_AR)

      testCase.assertTrue(QrcbMap.has_QRCID("QRCID1"))

      testCase.assertError(...
          @() QrcbMap.add("QRCID1", ~QRCB_AR), ...
          ?MException)

      QrcbMap.set("QRCID1", ~QRCB_AR)

      actQrcbAr = QrcbMap.get("QRCID1");
      testCase.assertEqual(actQrcbAr, ~QRCB_AR)
    end



    function test_add_false(T)
      QrcbMap = bicas.proc.QrcbMap(2);

      QrcbMap.add_false(string.empty(0, 1))

      QrcbMap.add_false(["QRCID_1"; "QRCID_2"])
      T.assertEqual(QrcbMap.get("QRCID_1"), false(2, 1))
      T.assertEqual(QrcbMap.get("QRCID_2"), false(2, 1))
    end



    function test_remove(testCase)
      QrcbMap = bicas.proc.QrcbMap(2);

      QrcbMap.add("QRCID1", logical([0; 1]))
      QrcbMap.add("QRCID2", logical([1; 0]))

      QrcbMap.remove("QRCID1")
      testCase.assertFalse(QrcbMap.has_QRCID("QRCID1"))
      testCase.assertTrue(QrcbMap.has_QRCID("QRCID2"))
    end



    function test_add_map___empty(testCase)
      QrcbMap      = bicas.proc.QrcbMap(3);
      AddedQrcbMap = bicas.proc.QrcbMap(3);
      ExpQrcbMap   = bicas.proc.QrcbMap(3);

      QrcbMap.add_map(AddedQrcbMap);

      testCase.assertEqual(QrcbMap, ExpQrcbMap)
    end



    function test_add_map___one_empty_1(testCase)
      QrcbMap      = bicas.proc.QrcbMap(3);
      QrcbMap.add(   "QRCID1", logical([0; 1; 0]));
      QrcbMap.add(   "QRCID2", logical([1; 0; 1]));

      AddedQrcbMap = bicas.proc.QrcbMap(3);

      ExpQrcbMap   = bicas.proc.QrcbMap(3);
      ExpQrcbMap.add("QRCID1", logical([0; 1; 0]));
      ExpQrcbMap.add("QRCID2", logical([1; 0; 1]));

      QrcbMap.add_map(AddedQrcbMap);

      testCase.assertEqual(QrcbMap, ExpQrcbMap)
    end



    function test_add_map___one_empty_2(testCase)
      QrcbMap      = bicas.proc.QrcbMap(3);

      AddedQrcbMap = bicas.proc.QrcbMap(3);
      AddedQrcbMap.add("QRCID1", logical([0; 1; 0]));
      AddedQrcbMap.add("QRCID2", logical([1; 0; 1]));

      ExpQrcbMap   = bicas.proc.QrcbMap(3);
      ExpQrcbMap.add(  "QRCID1", logical([0; 1; 0]));
      ExpQrcbMap.add(  "QRCID2", logical([1; 0; 1]));

      QrcbMap.add_map(AddedQrcbMap);

      testCase.assertEqual(QrcbMap, ExpQrcbMap)
    end



    function test_add_map___overlap(testCase)
      QrcbMap      = bicas.proc.QrcbMap(4);
      QrcbMap.add("QRCID1",              logical([0; 0; 0; 0]));
      QrcbMap.add("QRCID2_overlap",      logical([0; 0; 1; 1]));

      AddedQrcbMap = bicas.proc.QrcbMap(4);
      AddedQrcbMap.add("QRCID2_overlap", logical([0; 1; 0; 1]));
      AddedQrcbMap.add("QRCID3",         logical([1; 0; 1; 0]));

      ExpQrcbMap   = bicas.proc.QrcbMap(4);
      ExpQrcbMap.add("QRCID1",           logical([0; 0; 0; 0]));
      ExpQrcbMap.add("QRCID2_overlap",   logical([0; 1; 1; 1]));
      ExpQrcbMap.add("QRCID3",           logical([1; 0; 1; 0]));

      QrcbMap.add_map(AddedQrcbMap);

      testCase.assertEqual(QrcbMap, ExpQrcbMap)
    end



    function test_equality___empty(testCase)
      QrcbMap1a = bicas.proc.QrcbMap(1);
      QrcbMap1b = bicas.proc.QrcbMap(1);
      QrcbMap2  = bicas.proc.QrcbMap(2);

      testCase.assert_equal(    QrcbMap1a, QrcbMap1b)
      testCase.assert_not_equal(QrcbMap1a, QrcbMap2)
    end



    function test_equality___QRCB(testCase)
      QrcbMap1 = bicas.proc.QrcbMap(1);
      QrcbMap2 = bicas.proc.QrcbMap(1);

      QrcbMap1.add("QRCID1", logical([1]))
      testCase.assert_not_equal(QrcbMap1, QrcbMap2)

      QrcbMap2.add("QRCID1", logical([1]))
      testCase.assert_equal(    QrcbMap1, QrcbMap2)

      QrcbMap1.add("QRCID2", logical([2]))
      testCase.assert_not_equal(QrcbMap1, QrcbMap2)
    end



    function test_equality___unequal_QRCB(testCase)
      QrcbMap1 = bicas.proc.QrcbMap(2);
      QrcbMap2 = bicas.proc.QrcbMap(2);

      QrcbMap1.add("QRCID1", logical([0; 1]))
      QrcbMap2.add("QRCID1", logical([1; 1]))

      testCase.assert_not_equal(QrcbMap1, QrcbMap2)
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



      function assert_equal(testCase, QrcbMap1, QrcbMap2)
        testCase.assertTrue(isequal(QrcbMap1, QrcbMap2))
        testCase.assertEqual(       QrcbMap1, QrcbMap2)
        testCase.assertTrue(isequal(QrcbMap2, QrcbMap1))
        testCase.assertEqual(       QrcbMap2, QrcbMap1)
      end



      function assert_not_equal(testCase, QrcbMap1, QrcbMap2)
        testCase.assertFalse(isequal(QrcbMap1, QrcbMap2))
        testCase.assertNotEqual(     QrcbMap1, QrcbMap2)
        testCase.assertFalse(isequal(QrcbMap2, QrcbMap1))
        testCase.assertNotEqual(     QrcbMap2, QrcbMap1)
      end



  end    % methods(Access=private)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)
  end    % methods(Static, Access=private)



end
