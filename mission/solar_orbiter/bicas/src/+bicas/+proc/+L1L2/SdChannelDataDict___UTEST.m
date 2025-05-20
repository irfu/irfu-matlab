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
      SDID_ASR_AR = bicas.proc.L1L2.const.C.SDID_ASR_AR;

      SdcdDict = bicas.proc.L1L2.SdChannelDataDict();

      for i = 1:numel(SDID_ASR_AR)
        sdid = SDID_ASR_AR(i);

        Schd = bicas.proc.L1L2.SingleChannelData(...
          [NaN,2; 3,NaN; NaN,NaN]+i, ...
          logical([0; 1; 0]));
        SdcdDict.set(sdid, Schd);

        testCase.assertEqual(SdcdDict.nWholeRowIsNan, 1*i)
      end

      ExpSchd = bicas.proc.L1L2.SingleChannelData([NaN,2; 3,NaN; NaN,NaN]+3, logical([0; 1; 0]));
      ActSchd = SdcdDict.get(SDID_ASR_AR(3));

      testCase.assertEqual(ActSchd, ExpSchd)
      testCase.assertEqual(SdcdDict.nWholeRowIsNan, numel(SDID_ASR_AR)*1)
    end



  end    % methods(Test)



end
