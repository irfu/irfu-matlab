%
% UNFINISHED
%
% matlab.unittest automatic test code for bicas.proc.L1L2.ChannelDataAsrCollection.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ChannelDataAsrCollection___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_nRecords(testCase)
      Cdac = bicas.proc.L1L2.ChannelDataAsrCollection(0);
      testCase.assertEqual(Cdac.nRecords, 0)

      Cdac = bicas.proc.L1L2.ChannelDataAsrCollection(3);
      testCase.assertEqual(Cdac.nRecords, 3)
    end



    function test_set_get_nWholeRowIsNan(testCase)
      SDID_ASR_AR = bicas.proc.L1L2.const.C.SDID_ASR_AR;

      % 3x2
      SAMPLES_AR_0 = [NaN,2; 3,NaN; NaN,NaN];
      VSIB_AR      = logical([0; 1; 0]);

      Cdac = bicas.proc.L1L2.ChannelDataAsrCollection(3);

      for i = 1:numel(SDID_ASR_AR)
        sdid = SDID_ASR_AR(i);

        Schd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR_0+i, VSIB_AR);
        Cdac.set_channel(sdid, Schd);

        testCase.assertEqual(Cdac.nWholeRowIsNan, 1*i)
      end

      testCase.assertEqual(Cdac.nRecords, 3)

      ExpSchd = bicas.proc.L1L2.SingleChannelData(SAMPLES_AR_0+3, VSIB_AR);
      ActSchd = Cdac.get_channel(SDID_ASR_AR(3));

      testCase.assertEqual(ActSchd, ExpSchd)
      testCase.assertEqual(Cdac.nWholeRowIsNan, numel(SDID_ASR_AR)*1)
    end



  end    % methods(Test)



end
