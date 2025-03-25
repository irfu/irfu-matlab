%
% UNFINISHED
%
% matlab.unittest automatic test code for bicas.proc.L1L2.SdChannelDataDict.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SdChannelDataDict___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_set_get_nIsNan(testCase)
      SDID_ASR_AR = bicas.proc.L1L2.const.C.SDID_ASR_AR;

      SdcdDict = bicas.proc.L1L2.SdChannelDataDict();

      for i = 1:9
        sdid = SDID_ASR_AR(i);

        Sdcd = bicas.proc.L1L2.SdChannelData([NaN,2; 3,NaN]+i, logical([0; 1]));
        SdcdDict.set(sdid, Sdcd);

        testCase.assertEqual(SdcdDict.nIsNan, 2*i)
      end

      ExpSdcd = bicas.proc.L1L2.SdChannelData([NaN,2; 3,NaN]+3, logical([0; 1]));
      ActSdcd = SdcdDict.get(SDID_ASR_AR(3));

      testCase.assertEqual(ActSdcd, ExpSdcd)
      testCase.assertEqual(SdcdDict.nIsNan, 9*2)
    end



  end    % methods(Test)



end
