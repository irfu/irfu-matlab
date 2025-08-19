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
      Qrcbm = bicas.proc.QrcbMap(2);
      testCase.assertEqual(Qrcbm.qrcidAr, string.empty(0, 1))

      % IMPLEMENTATION NOTE: Deliberately add QRCIDs in non-alphanumeric order.
      Qrcbm.add("QRCID2", logical([0; 1]))
      testCase.assertEqual(Qrcbm.qrcidAr, ["QRCID2"])

      % Deliberately assert alphanumeric order for return value.
      Qrcbm.add("QRCID1", logical([1; 0]))
      testCase.assertEqual(Qrcbm.qrcidAr, ["QRCID1"; "QRCID2"])
    end



    function test_add_set_get_has_QRCID(testCase)
      QRCB_AR = logical([0; 1]);
      Qrcbm = bicas.proc.QrcbMap(2);

      testCase.assertError(...
        @() Qrcbm.set("QRCID1", ~QRCB_AR), ...
        ?MException)

      testCase.assertFalse(Qrcbm.has_QRCID("QRCID1"))

      Qrcbm.add("QRCID1", QRCB_AR)

      testCase.assertTrue(Qrcbm.has_QRCID("QRCID1"))

      testCase.assertError(...
        @() Qrcbm.add("QRCID1", ~QRCB_AR), ...
        ?MException)

      Qrcbm.set("QRCID1", ~QRCB_AR)

      actQrcbAr = Qrcbm.get("QRCID1");
      testCase.assertEqual(actQrcbAr, ~QRCB_AR)
    end



    function test_add_false(T)
      Qrcbm = bicas.proc.QrcbMap(2);

      Qrcbm.add_false(string.empty(0, 1))

      Qrcbm.add_false(["QRCID_1"; "QRCID_2"])
      T.assertEqual(Qrcbm.get("QRCID_1"), false(2, 1))
      T.assertEqual(Qrcbm.get("QRCID_2"), false(2, 1))
    end



    % function test_remove(testCase)
    %   Qrcbm = bicas.proc.QrcbMap(2);
    %
    %   Qrcbm.add("QRCID1", logical([0; 1]))
    %   Qrcbm.add("QRCID2", logical([1; 0]))
    %
    %   Qrcbm.remove("QRCID1")
    %   testCase.assertFalse(Qrcbm.has_QRCID("QRCID1"))
    %   testCase.assertTrue(Qrcbm.has_QRCID("QRCID2"))
    % end



    function test_union___empty(testCase)
      Qrcbm      = bicas.proc.QrcbMap(3);
      AddedQrcbm = bicas.proc.QrcbMap(3);
      ExpQrcbm   = bicas.proc.QrcbMap(3);

      Qrcbm.union(AddedQrcbm);

      testCase.assertEqual(Qrcbm, ExpQrcbm)
    end



    function test_union___one_empty_1(testCase)
      Qrcbm      = bicas.proc.QrcbMap(3);
      Qrcbm.add(   "QRCID1", logical([0; 1; 0]));
      Qrcbm.add(   "QRCID2", logical([1; 0; 1]));

      AddedQrcbm = bicas.proc.QrcbMap(3);

      ExpQrcbm   = bicas.proc.QrcbMap(3);
      ExpQrcbm.add("QRCID1", logical([0; 1; 0]));
      ExpQrcbm.add("QRCID2", logical([1; 0; 1]));

      Qrcbm.union(AddedQrcbm);

      testCase.assertEqual(Qrcbm, ExpQrcbm)
    end



    function test_union___one_empty_2(testCase)
      Qrcbm      = bicas.proc.QrcbMap(3);

      AddedQrcbm = bicas.proc.QrcbMap(3);
      AddedQrcbm.add("QRCID1", logical([0; 1; 0]));
      AddedQrcbm.add("QRCID2", logical([1; 0; 1]));

      ExpQrcbm   = bicas.proc.QrcbMap(3);
      ExpQrcbm.add(  "QRCID1", logical([0; 1; 0]));
      ExpQrcbm.add(  "QRCID2", logical([1; 0; 1]));

      Qrcbm.union(AddedQrcbm);

      testCase.assertEqual(Qrcbm, ExpQrcbm)
    end



    function test_union___overlap(testCase)
      Qrcbm      = bicas.proc.QrcbMap(4);
      Qrcbm.add("QRCID1",              logical([0; 0; 0; 0]));
      Qrcbm.add("QRCID2_overlap",      logical([0; 0; 1; 1]));

      AddedQrcbm = bicas.proc.QrcbMap(4);
      AddedQrcbm.add("QRCID2_overlap", logical([0; 1; 0; 1]));
      AddedQrcbm.add("QRCID3",         logical([1; 0; 1; 0]));

      ExpQrcbm   = bicas.proc.QrcbMap(4);
      ExpQrcbm.add("QRCID1",           logical([0; 0; 0; 0]));
      ExpQrcbm.add("QRCID2_overlap",   logical([0; 1; 1; 1]));
      ExpQrcbm.add("QRCID3",           logical([1; 0; 1; 0]));

      Qrcbm.union(AddedQrcbm);

      testCase.assertEqual(Qrcbm, ExpQrcbm)
    end



    function test_equality___empty(testCase)
      Qrcbm1a = bicas.proc.QrcbMap(1);
      Qrcbm1b = bicas.proc.QrcbMap(1);
      Qrcbm2  = bicas.proc.QrcbMap(2);

      testCase.assert_equal(    Qrcbm1a, Qrcbm1b)
      testCase.assert_not_equal(Qrcbm1a, Qrcbm2)
    end



    function test_equality___QRCB(testCase)
      Qrcbm1 = bicas.proc.QrcbMap(1);
      Qrcbm2 = bicas.proc.QrcbMap(1);

      Qrcbm1.add("QRCID1", logical([1]))
      testCase.assert_not_equal(Qrcbm1, Qrcbm2)

      Qrcbm2.add("QRCID1", logical([1]))
      testCase.assert_equal(    Qrcbm1, Qrcbm2)

      Qrcbm1.add("QRCID2", logical([2]))
      testCase.assert_not_equal(Qrcbm1, Qrcbm2)
    end



    function test_equality___unequal_QRCB(testCase)
      Qrcbm1 = bicas.proc.QrcbMap(2);
      Qrcbm2 = bicas.proc.QrcbMap(2);

      Qrcbm1.add("QRCID1", logical([0; 1]))
      Qrcbm2.add("QRCID1", logical([1; 1]))

      testCase.assert_not_equal(Qrcbm1, Qrcbm2)
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function assert_equal(testCase, Qrcbm1, Qrcbm2)
      testCase.assertTrue(isequal(Qrcbm1, Qrcbm2))
      testCase.assertEqual(       Qrcbm1, Qrcbm2)
      testCase.assertTrue(isequal(Qrcbm2, Qrcbm1))
      testCase.assertEqual(       Qrcbm2, Qrcbm1)
    end



    function assert_not_equal(testCase, Qrcbm1, Qrcbm2)
      testCase.assertFalse(isequal(Qrcbm1, Qrcbm2))
      testCase.assertNotEqual(     Qrcbm1, Qrcbm2)
      testCase.assertFalse(isequal(Qrcbm2, Qrcbm1))
      testCase.assertNotEqual(     Qrcbm2, Qrcbm1)
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
