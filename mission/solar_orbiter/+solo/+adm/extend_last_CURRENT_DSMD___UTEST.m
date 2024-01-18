%
% matlab.unittest automatic test code for
% solo.adm.extend_last_CURRENT_DSMD().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef extend_last_CURRENT_DSMD___UTEST < matlab.unittest.TestCase

  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test(DsmdArray1, timeExtensionDays, expDsmdArray)
        actDsmdArray = solo.adm.extend_last_CURRENT_DSMD(...
          DsmdArray1, timeExtensionDays);
        testCase.assertEqual(actDsmdArray, expDsmdArray)
      end

      function dt = dtu(varargin)
        dt = datetime(varargin{:}, 'TimeZone', 'UTCLeapSeconds');
      end

      EDV         = solo.adm.DSMD.empty(0,1);
      CURRENT_DSI = 'SOLO_L1_RPW-BIA-CURRENT';
      CWF_DSI     = 'SOLO_L2_RPW-LFR-SURV-CWF-E';

      CUR1a = solo.adm.DSMD('path', CURRENT_DSI, 99, true, dtu(2000, 1, 1, 0, 0, 0), dtu(2001, 2, 02, 12, 34, 56));
      CUR1b = solo.adm.DSMD('path', CURRENT_DSI, 99, true, dtu(2000, 1, 1, 0, 0, 0), dtu(2001, 2, 16, 12, 34, 56));
      CUR2a = solo.adm.DSMD('path', CURRENT_DSI, 99, true, dtu(2002, 1, 1, 0, 0, 0), dtu(2002, 2, 02, 12, 34, 56));
      CUR2b = solo.adm.DSMD('path', CURRENT_DSI, 99, true, dtu(2002, 1, 1, 0, 0, 0), dtu(2002, 2, 16, 12, 34, 56));

      CWF1  = solo.adm.DSMD('path', CWF_DSI, 99, true, dtu(2000, 2, 3, 0, 0, 0), dtu(2001, 2, 02, 12, 34, 56));
      CWF2  = solo.adm.DSMD('path', CWF_DSI, 99, true, dtu(2001, 2, 3, 0, 0, 0), dtu(2002, 2, 16, 12, 34, 56));

      test(EDV,   14, EDV);
      test(CUR1a, 14, CUR1b);
      test(CWF1,  14, CWF1);
      test([CUR1a; CWF1], 14, [CUR1b; CWF1]);
      test([CUR1a; CWF1; CUR2a; CWF2], 14, [CUR1a; CWF1; CUR2b; CWF2]);
      test([CUR2a; CWF1; CUR1a; CWF2], 14, [CUR2b; CWF1; CUR1a; CWF2]);
    end



  end    % methods(Test)

end
