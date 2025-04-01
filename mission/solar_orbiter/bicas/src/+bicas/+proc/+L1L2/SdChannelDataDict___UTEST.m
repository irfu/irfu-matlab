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



    function test_set_get_nWholeRowIsNan(testCase)
      SDID_ASR_AR = bicas.proc.L1L2.const.C.SDID_ASID_DICT.keys;

      SdcdDict = bicas.proc.L1L2.SdChannelDataDict();

      for i = 1:numel(SDID_ASR_AR)
        sdid = SDID_ASR_AR(i);

        Sdcd = bicas.proc.L1L2.SdChannelData(...
          [NaN,2; 3,NaN; NaN,NaN]+i, ...
          logical([0; 1; 0]));
        SdcdDict.set(sdid, Sdcd);

        testCase.assertEqual(SdcdDict.nWholeRowIsNan, 1*i)
      end

      ExpSdcd = bicas.proc.L1L2.SdChannelData([NaN,2; 3,NaN; NaN,NaN]+3, logical([0; 1; 0]));
      ActSdcd = SdcdDict.get(SDID_ASR_AR(3));

      testCase.assertEqual(ActSdcd, ExpSdcd)
      testCase.assertEqual(SdcdDict.nWholeRowIsNan, numel(SDID_ASR_AR)*1)
    end



  end    % methods(Test)



end
