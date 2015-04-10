classdef test_mms_spinfit < matlab.unittest.TestCase
  %TEST_MMS_DEFATT_PHASE
  
  properties
  end
  
  methods (Test)
    function test_constant_phase(testCase)
      %% Test constant spin
      spinRate = 7.3; %RPM
      Ex = 0.3; Ey = 0.5; adcOff_12 = -0.45; adcOff_34 = 0.65;
      t0 = int64(481744867369743068); % 2015-04-08T06:00:00.185743068Z
      phiDeg0 = 267.5688;
      MMS_CONST=mms_constants;
      
      %0.5 sec sampling DEFATT
      timeSecPhase = (1:3600)'*.5;
      defatt.time = int64(timeSecPhase*1e9) + t0;
      defatt.zphase = mod(phiDeg0 + 360*timeSecPhase*spinRate/60,360);
      %32 sps sampling E
      timeSecReq = (1:(1800*32))'/32;
      timeEpochTT200Req = int64(timeSecReq*1e9) + t0;
      phaseComputed = mms_defatt_phase(defatt,timeEpochTT200Req);
      phaseRad = unwrap(phaseComputed.data*pi/180);
      phaseRad12 = phaseRad - MMS_CONST.Phaseshift.e12;
      phaseRad34 = phaseRad - MMS_CONST.Phaseshift.e34;
      e12 = Ex*cos(phaseRad12) + Ey*sin(phaseRad12) + adcOff_12;
      e34 = Ex*cos(phaseRad34) + Ey*sin(phaseRad34) + adcOff_34;
      
      % Plot
      %irf_plot([EpochTT2000(timeEpochTT200Req).toEpochUnix().epoch e12 e34])
      
      % Check with despin
      %dE = mms_sdp_despin(e12-adcOff_12,e34-adcOff_34,phaseComputed.data);
      %irf_plot([EpochTT2000(timeEpochTT200Req).toEpochUnix().epoch dE])
      %legend('Ex','Ey')
      
      % Compute spinfit for E12
      [~,sfit,~,~,~] = ...
                c_efw_spinfit_mx(3,10,5,...
                double(timeEpochTT200Req-t0)'*1e-9,e12',phaseRad12');
      sfit(sfit==-1.5900e+09)=NaN; sfit = sfit'; sfit(isnan(sfit(:,1)),:)=[];
      
      %Test 12
      testCase.verifyEqual(median(sfit(:,1)),adcOff_12,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,2)),Ex,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,3)),Ey,'AbsTol',1e-6);
      
      % Compute spinfit for E34
      [~,sfit,~,~,~] = ...
                c_efw_spinfit_mx(3,10,5,...
                double(timeEpochTT200Req-t0)'*1e-9,e34',phaseRad34');
      sfit(sfit==-1.5900e+09)=NaN; sfit = sfit'; sfit(isnan(sfit(:,1)),:)=[];
      
      %Test 34
      testCase.verifyEqual(median(sfit(:,1)),adcOff_34,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,2)),Ex,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,3)),Ey,'AbsTol',1e-6);
    end
  end
end