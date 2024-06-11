function [NeScp, NeScpQualityBit, codeVerStr] = psp2ne(PSP)
%SOLO.PSP2NE  Convert probe-to-spacecraft potential to electron density
%
% [NeScp, NeScpQualityBit, codeVerStr] = solo.psp2ne(PSP)
%
% Convert probe-to-spacecraft (PSP) potential to electron density (NeScp).
%
% The calibration is based on the RPW/QTN/Fpe data.
%
%
% RETURN VALUES
% =============
% NeScp
%       Electron density (derived from "SCP", hence the name).
% NeScpQualityBit
%       Binary value that specifies whether the density value
%       seems bad or not. 1=Bad, 0=Can not find any problem.
%       Must not be NaN. (Currently (2023-08-10) not sure if this
%       is strictly in agreement with conventions, but that is
%       what BICAS requires).
% codeVerStr
%       Date-formatted version string for the code that implements the function
%       *code*, including calibration data. This is used by BICAS for setting
%       the relevant CDF global attribute to indicate the version of the
%       algorithm used to produce a particular dataset (CDF file).
%       Must be on the form of a human-readable UTC timestamp string which
%       conforms to regular expression bicas.proc.L2L3.ext.CODE_VER_STR_REGEXP.
%
% NOTE: This function is used by BICAS for producing official L3 datasets. It
%       must therefore have an interface (name, arguments, return values) that
%       is compatible with BICAS.
%
% Calibration using plasma line
% see Dropbox/Solar_Orbiter/Science data/InFlight Cal/Ncalpsp2ne_calibrate.m

Cal = [];

%===============================================================================
% Timestamp string that represent the version of the function. This string is
% used by BICAS to set a CDF global attribute in official datasets for
% traceability.
%
% IMPORTANT NOTE: This value is meant to be be updated manually to
% the approximate current date when the calibration data or algorithm is
% updated (or possibly when calibration data or algorithm risk being
% unintentionally changed due to refactoring).
% * It should *NOT* be an automatically set timestamp (e.g. current time).
% * It should *NOT* be updated for unrelated code changes, e.g. comments or
%   variable name changes.
%===============================================================================
codeVerStr = '2023-12-19T15:11:00';

AddEntry('2020-03-08T00:00:00Z/2020-05-18T04:05:54Z',[0.8889  3.4389]); % Based on data from 2020-04-07
AddEntry('2020-05-18T04:05:55Z/2020-05-29T23:59:59Z',[0.8154  4.5562]);
AddEntry('2020-05-30T00:00:00Z/2020-08-11T21:27:02Z',[0.5310  4.4010]); % Based up to 2020-07-05T23:59:59Z
AddEntry('2020-08-11T21:27:03Z/2020-09-04T23:59:59Z',[0.6593  3.8785]); % Data until Sept 4
AddEntry('2020-09-05T00:00:00Z/2020-09-13T23:59:59Z',[0.7148  3.6008]); % Data for Sept 5
AddEntry('2020-09-14T00:00:00Z/2020-11-02T23:59:59Z',[0.6726  3.3577]);
AddEntry('2020-11-03T00:00:00Z/2020-11-30T23:59:59Z',[0.7367  3.5170]);
AddEntry('2020-12-01T00:00:00Z/2020-12-14T23:59:59Z',[0.6834  3.8871]);
AddEntry('2020-12-15T00:00:00Z/2020-12-31T23:59:59Z',[0.5556  3.8249]);
%%%%%%%%%%%%%%%%%%%%%%%%%2021%%%%%%%%%%%%%%%%%%%%%%%%%%
AddEntry('2021-01-01T00:00:00Z/2021-01-06T23:59:59Z',[0.5905  4.0923]);
AddEntry('2021-01-07T00:00:00Z/2021-01-12T05:49:59Z',[0.6730  4.1837]); %2
AddEntry('2021-01-12T05:50:00Z/2021-01-17T23:59:59Z',[0.7462  4.5630]); %3
AddEntry('2021-01-18T00:00:00Z/2021-02-16T23:59:59Z',[0.3524  3.8309]); %4
AddEntry('2021-02-17T00:00:00Z/2021-03-20T01:29:59Z',[0.7651  3.4971]); %5
AddEntry('2021-03-20T01:30:00Z/2021-03-22T19:29:59Z',... %6
  [0.6460 + 3.5047i   0.3899 + 3.7683i],1.0297);
AddEntry('2021-03-22T19:30:00Z/2021-04-04T03:59:59Z',[0.7884  3.3714]); %7
AddEntry('2021-04-04T04:00:00Z/2021-07-25T23:59:59Z',...%8
  [0.7125 + 3.0114i   0.4926 +  3.2371i], 1.050);
AddEntry('2021-07-26T00:00:00Z/2021-08-02T23:59:59Z',[0.7396  3.2209]); %9
AddEntry('2021-08-03T00:00:00Z/2021-08-05T07:59:59Z',[0.7694  3.3844]); %10
AddEntry('2021-08-05T08:00:00Z/2021-08-08T23:59:59Z',[0.6615  3.1782]); %11
AddEntry('2021-08-09T00:00:00Z/2021-08-11T23:59:59Z',[0.6802  3.2830]); %12
AddEntry('2021-08-12T00:00:00Z/2021-08-29T19:59:59Z',... %13
  [0.7572 + 3.3979i  0.4088 + 3.6251i],0.6522);
AddEntry('2021-08-29T18:00:00Z/2021-08-31T23:59:59Z',... %14
  [0.6898 + 3.7050i  0.4299 + 3.6844i],-0.0793);
AddEntry('2021-09-01T00:00:00Z/2021-09-04T00:59:59Z',... %15
  [0.6245 + 3.5472i  0.2589 + 3.8734i],0.8925);
AddEntry('2021-09-04T01:00:00Z/2021-09-11T23:59:59Z',[0.6916  3.5999]); %16
AddEntry('2021-09-12T00:00:00Z/2021-09-27T22:05:45Z',... %17
  [0.5830 + 3.6131i  0.4904 + 3.5927i],0.6390);
AddEntry('2021-09-27T22:05:46Z/2021-10-05T01:59:59Z',... %18
  [0.7342 + 3.5247i  0.4475 + 3.4873i],-0.1302);
AddEntry('2021-10-05T02:00:00Z/2021-10-09T23:59:59Z',[0.7239  3.3958]); %19
AddEntry('2021-10-10T00:00:00Z/2021-10-14T23:59:59Z',[0.8061  3.2660]); %20
AddEntry('2021-10-15T00:00:00Z/2021-10-15T23:59:59Z',[0.6726  3.2377]); %21
AddEntry('2021-10-16T00:00:00Z/2021-10-23T23:59:59Z',... %22
  [0.8444 + 3.1733i  0.5509 + 3.3010i], 0.550);
AddEntry('2021-10-24T00:00:00Z/2021-11-05T23:59:59Z',... %23
  [0.6579 + 3.1300i  0.4599 + 3.2233i], 0.4709);
AddEntry('2021-11-06T00:00:00Z/2021-11-20T23:59:59Z',[0.7547  2.8010]); %24
AddEntry('2021-11-21T00:00:00Z/2021-12-18T23:59:59Z',... %25
  [0.6849 + 2.7147i  0.4785 + 2.9694i], 1.1);
AddEntry('2021-12-19T00:00:00Z/2021-12-31T23:59:59Z',[0.5435  2.8046]); %26
%%%%%%%%%%%%%%%%%%%%%%%%%2022%%%%%%%%%%%%%%%%%%%%%%%%%%


AddEntry('2022-01-01T00:00:00Z/2022-01-09T23:59:59Z',... %27
  [0.8065 + 2.7517i  0.4831 + 2.8837i],0.4081);
AddEntry('2022-01-10T00:00:00Z/2022-01-31T23:59:59Z',[0.6289  2.9140]); %28
AddEntry('2022-02-01T00:00:00Z/2022-02-07T23:59:59Z',[0.7092  3.0180]); %29
AddEntry('2022-02-08T00:00:00Z/2022-02-14T23:59:59Z',[0.5814  3.2496]); %30
AddEntry('2022-02-15T00:00:00Z/2022-02-16T23:59:59Z',... %31
  [0.6536 + 3.1604i  0.3484 + 3.4695i],1.3932);
AddEntry('2022-02-17T00:00:00Z/2022-02-21T14:29:59Z',[0.4545  3.3593]); %32
AddEntry('2022-02-21T14:30:00Z/2022-02-21T23:59:59Z',[0.8621  3.6587]); %33
AddEntry('2022-02-22T00:00:00Z/2022-02-24T13:29:59Z',[0.6452  3.1224]); %34
AddEntry('2022-02-24T13:30:00Z/2022-03-01T23:59:59Z',... %35
  [0.4975 + 3.3311i 0.1812 + 3.8343i],1.5912);
AddEntry('2022-03-02T00:00:00Z/2022-03-02T23:59:59Z',[0.8000  3.6187]); %36
AddEntry('2022-03-03T00:00:00Z/2022-03-07T23:59:59Z',[0.6579  3.9701]); %37
AddEntry('2022-03-08T00:00:00Z/2022-03-08T05:29:59Z',[0.3676  3.7905]); %38
AddEntry('2022-03-08T05:30:00Z/2022-03-08T14:45:59Z',... %39
  [0.5917 + 3.6481i  0.3344 + 3.8349i],0.7302);
AddEntry('2022-03-08T14:46:00Z/2022-03-10T21:59:59Z',... %40
  [0.5682 + 3.9543i  0.0579 + 4.5401i],1.1516);
AddEntry('2022-03-10T22:00:00Z/2022-03-12T07:59:59Z',[0.4854  3.9150]); %41
AddEntry('2022-03-12T08:00:00Z/2022-03-13T06:59:59Z',[0.4902  3.9480]); %42
AddEntry('2022-03-13T07:00:00Z/2022-03-15T23:59:59Z',[0.5682  4.2383]); %43
AddEntry('2022-03-16T00:00:00Z/2022-03-30T14:59:59Z',[0.4695  4.3045]); %44
AddEntry('2022-03-30T15:00:00Z/2022-04-04T23:59:59Z',[0.5155  3.7450]); %45
AddEntry('2022-04-05T00:00:00Z/2022-04-11T23:59:59Z',[0.2427 3.6602]); %46
AddEntry('2022-04-12T00:00:00Z/2022-04-20T23:59:59Z',... %47
  [0.3460 + 3.0611i  0.1776 + 3.5813i],3.0865);
AddEntry('2022-04-21T00:00:00Z/2022-05-15T23:59:59Z',[0.6803 2.0541]); %48
AddEntry('2022-05-16T00:00:00Z/2022-05-28T23:59:59Z',... %49
  [0.7353 + 1.6034i  0.2667 + 2.6383i],2.2049);
AddEntry('2022-05-29T00:00:00Z/2022-08-08T23:59:59Z',... %50
  [0.9709 + 0.8459i  0.4202 + 2.0656i],2.2271);
AddEntry('2022-08-09T00:00:00Z/2022-09-06T09:59:59Z',... %51
  [0.5780 + 1.9516i  0.2519 + 2.7794i],2.5448);
AddEntry('2022-09-06T10:00:00Z/2022-09-07T07:59:59Z',... %52
  [0.8065 + 1.6715i  0.5263 + 2.2597i],2.0973);
AddEntry('2022-09-07T08:00:00Z/2022-09-14T06:59:59Z',[0.6757 2.2607]); %53
AddEntry('2022-09-14T07:00:00Z/2022-10-18T19:59:59Z',... %54
  [0.8547 + 2.9232i  0.2128 + 3.8295i],1.7300);
AddEntry('2022-10-18T20:00:00Z/2022-11-03T23:59:59Z',[0.3968 3.1672]); %55
AddEntry('2022-11-04T00:00:00Z/2022-11-29T23:59:59Z',... %56
  [0.4405 + 2.4799i  0.2083 + 3.1493i],3.2500);
AddEntry('2022-11-30T00:00:00Z/2022-12-09T18:15:59Z',... %57
  [0.6897 + 1.7440i  0.2457 + 2.9565i],2.7252);
AddEntry('2022-12-09T18:16:00Z/2022-12-19T01:04:59Z',... %58
  [0.5208 + 1.9359i  0.2242 + 2.7081i],2.6131);
AddEntry('2022-12-19T01:05:00Z/2022-12-20T15:04:39Z',[0.4545 1.9587]); %59
AddEntry('2022-12-20T15:04:40Z/2022-12-31T23:59:59Z',[0.7194 1.5239]); %60

%% calibrate
CalR = Cal.resample(PSP);
NeScp = PSP;


NeScp.data = exp(CalR.x.data.*NeScp.data + CalR.y.data);


timeOutsideInterval = irf_time('2022-12-31T23:39:59Z','utc>ttns');
NeScp.data(NeScp.time.epoch > timeOutsideInterval)= NaN;



NeScp.name = 'NeScp';
NeScp.units = 'cm^-3';
NeScp.siConversion = 'cm^-3>1e6*m^-3';
NeScp.userData = '';

% NOTE: Setting temporary (but legal) return value for return variable that is
% not yet used by BICAS (2023-08-10).
% NOTE: Overwrite every value with zero in order to also overwrite Nan which
% may otherwise inherited from NeScp.
% NOTE: Density from TNR plasma line used to calibrate NeScp only measures up
% to 122 cc, everything above that value is uncertain, therefore is flagged.
% Low values of NeScp i.e <2 cc are also uncertain.
NeScpQualityBit = TSeries(NeScp.time, ones(size(NeScp.data)));
NeScpQualityBit.data(NeScp.data<=122 & NeScp.data>=2) = 0;
NeScpQualityBit.data(isnan(NeScpQualityBit.data)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Help function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function AddEntry(TintS, calData, PSPintersection)
    %Add new calibration entry
    % For two-fit calibration use three input arguments

    if ~isreal(calData(1)) && (nargin<3 || isempty(PSPintersection))
      errS = ['Invalid two-fit cal entry at: ' TintS];
      irf.log('critical',errS)
      error(errS)
    end

    CalEntry = irf.ts_vec_xy(irf.tint(TintS),repmat(calData,2,1));

    if nargin>2 && ~isempty(PSPintersection)

      checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

      if ~isempty(checkInterval)
        [CalEntry] = TwoFitCalibration(checkInterval,PSPintersection,CalEntry);
      else
        CalEntry.data(1:end,2) = imag(CalEntry.x.data);
        CalEntry.data(1:end,1) = real(CalEntry.x.data);
      end
    end

    if isempty(Cal), Cal = CalEntry;
    else, Cal = Cal.combine(CalEntry);
    end

    function [C] = TwoFitCalibration(PSPint,y_eq,CalData)
      CalData = CalData.resample(PSPint);

      %Identify the data points corresponding to each fit
      CalR1 = CalData(PSPint.data<y_eq).x;
      CalR2 = CalData(PSPint.data>=y_eq).y;
      %Identify NaNs
      CalRnan= CalData(isnan(PSPint.data));

      %Create a TSeries with the data points of each fit
      C1 = irf.ts_vec_xy(CalR1.time,[real(CalR1.data) imag(CalR1.data)]);
      C2 = irf.ts_vec_xy(CalR2.time,[real(CalR2.data) imag(CalR2.data)]);
      Cnan = irf.ts_vec_xy(CalRnan.time,[ones(length(CalRnan.data),1) ones(length(CalRnan.data),1)]);


      %Merge both fits and NaNs to keep the same length as the input
      %If statements added for robustness. An error will be send if
      %combining empty objects.
      if ~isempty(C1) && ~isempty(C2)
        C = C1.combine(C2);
      elseif ~isempty(C1) && isempty(C2)
        C = C1;
      elseif ~isempty(C2) && isempty(C1)
        C = C2;
      elseif isempty(C1) && isempty(C2) && ~isempty(Cnan)
        C = Cnan;
      elseif isempty(C1) && isempty(C2) && isempty(Cnan)
        C = CalData;
        irf.log('critical', 'no data at all ?!?')
      end

      if ~isempty(Cnan)
        C = C.combine(Cnan);
      end


    end
  end
end
