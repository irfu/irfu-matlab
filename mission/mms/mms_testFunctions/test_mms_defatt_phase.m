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
    function test_spinup(testCase)
      % Add an increase of spinRate by FACTOR in the end of interval
      spinRate = 3.1; %RPM
      factor = 1.1;
      %0.5 sec sampling DEFATT
      timeSecPhase = (1:3600)'*.5;
      defatt.time = int64(timeSecPhase*1e9);
      defatt.zphase = 360*timeSecPhase*spinRate/60;
      tSpinup = timeSecPhase(end-120); p0 = defatt.zphase(end-120);
      iiZ = timeSecPhase>=tSpinup;
      defatt.zphase(iiZ) =  p0 + 360*(timeSecPhase(iiZ)-tSpinup)*factor*spinRate/60;
      defatt.zphase = mod(defatt.zphase,360);
      %32 sps sampling E
      timeSecReq = (1:(1800*32))'/32;
      timeEpochTT200Req = int64(timeSecReq*1e9);
      phaseExpected = 360*timeSecReq*spinRate/60;
      iiE = timeSecReq>=tSpinup;
      phaseExpected(iiE) =  p0 + 360*(timeSecReq(iiE)-tSpinup)*factor*spinRate/60;
      phaseExpected = mod(phaseExpected,360);
      phaseComputed= mms_defatt_phase(defatt,timeEpochTT200Req);
      % Unwrap to continous phase to avoid issues with 0 ~= 359.9999... deg.
      testCase.verifyEqual(unwrap(pi/180*(phaseComputed.data)),unwrap(pi/180*(phaseExpected)),'AbsTol',pi/180*(1e-6));
      diffangle = phaseComputed.data - phaseExpected;
      diffangle = min([diffangle';360-diffangle']);
      testCase.verifyEqual(diffangle,diffangle*0,'AbsTol',1e-6);
    end
  end
end