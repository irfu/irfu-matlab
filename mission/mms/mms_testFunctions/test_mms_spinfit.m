classdef test_mms_spinfit < matlab.unittest.TestCase
  %TEST_MMS_DEFATT_PHASE

  properties
  end

  methods (Test)
    function test_c_efw_spinfit_mx_simple(testCase)
      %% Test constant spin
      spinRate = 7.3; %RPM
      Ex = 0.3; Ey = 0.5; adcOff_12 = -0.45; adcOff_34 = 0.65;
      t0 = int64(481744867369743068); % 2015-04-08T06:00:00.185743068Z
      phiDeg0 = 267.5688;
      MMS_CONST=mms_constants;
      DEBUG = false; % Show some plots?
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

      % Check with despin
      dE = mms_sdp_despin(e12-adcOff_12,e34-adcOff_34,phaseComputed.data); %#ok<NASGU>

      if(DEBUG)
        % Plot
        irf_plot([EpochTT(timeEpochTT200Req).epochUnix e12 e34]); %#ok<UNRCH>
        irf_plot([EpochTT(timeEpochTT200Req).epochUnix dE])
        legend('Ex','Ey')
      end

      % Compute spinfit for E12
      [~,sfit,~,~,~] = ...
        c_efw_spinfit_mx(3,10,3,...
        double(timeEpochTT200Req-t0)'*1e-9,e12',phaseRad12');
      sfit(sfit==-1.5900e+09)=NaN; sfit = sfit'; sfit(isnan(sfit(:,1)),:)=[];

      %Test 12
      testCase.verifyEqual(median(sfit(:,1)),adcOff_12,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,2)),Ex,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,3)),Ey,'AbsTol',1e-6);

      % Compute spinfit for E34
      [~,sfit,~,~,~] = ...
        c_efw_spinfit_mx(3,10,3,...
        double(timeEpochTT200Req-t0)'*1e-9,e34',phaseRad34');
      sfit(sfit==-1.5900e+09)=NaN; sfit = sfit'; sfit(isnan(sfit(:,1)),:)=[];

      %Test 34
      testCase.verifyEqual(median(sfit(:,1)),adcOff_34,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,2)),Ex,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,3)),Ey,'AbsTol',1e-6);
    end
    function test_mms_spinfit_mx_sinus_with_offset(testCase)
      % Simplest of test cases. Input data is only a pure sinus wave with
      % one static offset, a static spinrate is also applied. No
      % phaseshift or anything complicated.
      DEBUG = false; % Show some plots...
      fitEvery = 5e9; fitInterv = 20e9; %ns
      maxIt = 3; minPts = 10; nTerms = 3;
      spinRate = 3.1; %rpm
      % 0.5 sec sampling
      timeSecPhase = (0:100)'*.5;
      timeTT2000 = (timeSecPhase*1e9); %ns
      t0 = 5e9; % first even fitEvery, accounting for leap seconds...
      zphase = mod(360*timeSecPhase*spinRate/60, 360); % deg.
      radPhase = unwrap(zphase*pi/180); % rad, unwrapped
      % Pure sinus amplitude 1 with offset + 2
      expected_A = 2;
      expected_B = 0;
      expected_C = 1;
      % Data = A + B*cos(wt) + C*sin(wt)
      dataInput = expected_A + expected_B*cos(radPhase) + expected_C*sin(radPhase);
      % Call mms_spinfit_m, which in turn call the mex function.
      [~, sfit, ~, ~, ~] = mms_spinfit_m( maxIt, minPts, nTerms, ...
        timeTT2000, dataInput, radPhase, fitEvery, fitInterv, t0 );
      % Remove bad values
      sfit(isnan(sfit(:,1)),:) = [];
      diffA = median(sfit(:,1)) - expected_A;
      diffB = median(sfit(:,2)) - expected_B;
      diffC = median(sfit(:,3)) - expected_C;
      if(DEBUG)
        timeEpoch = irf_time(int64(timeTT2000),'ttns>epoch'); %#ok<UNRCH>
        extrapResult = median(sfit(:,1)) + median(sfit(:,2))*cosd(zphase) + median(sfit(:,3))*sind(zphase);
        figure; h=irf_plot({[timeEpoch, zphase],[timeEpoch, dataInput],[timeEpoch, extrapResult]});
        title(h(1),'Comparison of phase, input data and extrapolated data from median of spinfit coeff.');
        ylabel(h(1),{'Zphase','[deg]'})
        ylabel(h(2),{'Input data','y=A + ...','B*cosd(zphase) + ...',...
          'C*sind(zphase)','[mV/m]'})
        ylabel(h(3),{'Extrapol data from median spinfit',...
          'y=sfit(1) + ...','sfit(2)*cosd(zphase) + ...',...
          'sfit)*sind(zphase)','[mV/m]'}) % actually "median(sfit(1,:))" but write something shorter
      end
      % Verify all fits are as expected
      testCase.verifyEqual(diffA, diffA*0,'AbsTol',1e-6);
      testCase.verifyEqual(diffB, diffB*0,'AbsTol',1e-6);
      testCase.verifyEqual(diffC, diffC*0,'AbsTol',1e-6);
    end
    function test_mms_spinfit_mx_simple(testCase)
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
      %irf_plot([EpochTT(timeEpochTT200Req).epochUnix e12 e34])

      % Check with despin
      %dE = mms_sdp_despin(e12-adcOff_12,e34-adcOff_34,phaseComputed.data);
      %irf_plot([EpochTT(timeEpochTT200Req).epochUnix dE])
      %legend('Ex','Ey')

      % Compute spinfit for E12
      fitEvery = 5e9; fitInterv = 20e9; %ns
      [~, sfit, ~, ~, ~] = mms_spinfit_m( 3, 10, 5, ...
        timeEpochTT200Req, e12, phaseRad12, fitEvery, fitInterv, t0 );

      sfit(sfit==-1.5900e+09)=NaN; sfit(isnan(sfit(:,1)),:)=[];

      %Test 12
      testCase.verifyEqual(median(sfit(:,1)),adcOff_12,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,2)),Ex,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,3)),Ey,'AbsTol',1e-6);

      % Compute spinfit for E34
      [~, sfit, ~, ~, ~] = mms_spinfit_m( 3, 10, 5, ...
        timeEpochTT200Req, e34, phaseRad34, fitEvery, fitInterv, t0 );

      sfit(sfit==-1.5900e+09)=NaN; sfit(isnan(sfit(:,1)),:)=[];

      %Test 34
      testCase.verifyEqual(median(sfit(:,1)),adcOff_34,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,2)),Ex,'AbsTol',1e-6);
      testCase.verifyEqual(median(sfit(:,3)),Ey,'AbsTol',1e-6);
    end
  end
end
