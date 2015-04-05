classdef test_mms_defatt_phase < matlab.unittest.TestCase
  %TEST_MMS_DEFATT_PHASE
  
  properties
  end
  
  methods (Test)
    function test_constant_phase(testCase)
      spinRate = 3.1; %RPM
      %0.5 sec sampling DEFATT
      timeSecPhase = (1:3600)'*.5;
      defatt.time = int64(timeSecPhase*1e9);
      defatt.zphase = mod(360*timeSecPhase*spinRate/60,360);
      %32 sps sampling E
      timeSecReq = (1:(1800*32))'/32;
      timeEpochTT200Req = int64(timeSecReq*1e9);
      phaseExpected = mod(360*timeSecReq*spinRate/60,360);
      phaseComputed= mms_defatt_phase(defatt,timeEpochTT200Req);
      diffangle = phaseComputed.data - phaseExpected;
      diffangle = min([diffangle';360-diffangle']);
      testCase.verifyEqual(diffangle,diffangle*0,'AbsTol',1e-6);
    end
  end
end