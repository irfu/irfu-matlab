function [NeScp, codeVerStr] = psp2ne(PSP)
%SOLO.PSP2NE  Convert probe-to-spacecraft potential to electron density
%
% [NeScp,codeVerStr] = solo.psp2ne(PSP)
%
% Convert probe-to-spacecraft (PSP) potential to electron density (NeScp)
%
% The calibration is based on the RPW/QTN/Fpe data
%
% Outputs:
%   NeScp      - Electron density
%   codeVerStr - Version string. Used by BICAS.
%
% NOTE: This function is used by BICAS for producing official datasets.
%
% Calibration using plasma line 
% see Dropbox/Solar_Orbiter/Science data/InFlight Cal/Ncalpsp2ne_calibrate.m

Cal = [];

%===========================================================================
% Date string that represent the version of the function. This string is
% used by BICAS to set a CDF global attribute in official datasets for
% traceability.
% NOTE: This value is meant to be be updated by hand, not by an automatic
% timestamp, so that a constant value represents the same algorithm.
%===========================================================================
codeVerStr = '2023-07-17T11:38:00';

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

%======================2 fits=========================%
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
AddEntry('2022-03-12T08:00:00Z/2022-03-12T23:59:59Z',[0.5102  3.9640]); %42
AddEntry('2022-03-13T00:00:00Z/2022-03-31T23:59:59Z',[0.4049  4.0706]); %43
AddEntry('2022-04-12T00:00:00Z/2022-04-20T23:59:59Z',... %4~
  [0.3460 + 3.0611i  0.1776 + 3.5813i],3.0865); 

%% calibrate
CalR = Cal.resample(PSP);
NeScp = PSP; 


NeScp.data = exp(CalR.x.data.*NeScp.data + CalR.y.data);


timeOutsideInterval = irf_time('2022-04-20T23:59:59Z','utc>ttns');
NeScp.data(NeScp.time.epoch > timeOutsideInterval)= NaN;
    


NeScp.name = 'NeScp';
NeScp.units = 'cm^-3';
NeScp.siConversion = 'cm^-3>1e6*m^-3';
NeScp.userData = '';

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

