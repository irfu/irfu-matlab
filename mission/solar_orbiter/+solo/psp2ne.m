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



%===========================================================================
% Date string that represent the version of the function. This string is
% used by BICAS to set a CDF global attribute in official datasets for
% traceability.
% NOTE: This value is meant to be be updated by hand, not by an automatic
% timestamp, so that a constant value represents the same algorithm.
%===========================================================================
codeVerStr = '2022-05-27T16:04:00';


% Based on data from 2020-04-07
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-03-08T00:00:00Z/2020-05-18T04:05:54Z'),...
  repmat([0.8889   3.4389],2,1));
Cal = CalEntry;

CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-05-18T04:05:55Z/2020-05-29T23:59:59Z'),...
  repmat([0.8154   4.5562],2,1));
Cal = Cal.combine(CalEntry);

% cal based up to 2020-07-05T23:59:59Z
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-05-30T00:00:00Z/2020-08-11T21:27:02Z'),...
  repmat([0.5310   4.4010],2,1));
Cal = Cal.combine(CalEntry);

% data until Sept 4
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-08-11T21:27:03Z/2020-09-04T23:59:59Z'),...
  repmat([0.6593  3.8785],2,1));
Cal = Cal.combine(CalEntry);

% data for Sept 5
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-09-05T00:00:00Z/2020-09-13T23:59:59Z'),...
  repmat([0.7148  3.6008],2,1));
Cal = Cal.combine(CalEntry);

%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-09-14T00:00:00Z/2020-11-02T23:59:59Z'),...
  repmat([0.6726  3.3577],2,1));
Cal = Cal.combine(CalEntry);

%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-11-03T00:00:00Z/2020-11-30T23:59:59Z'),...
  repmat([0.7367  3.5170],2,1));
Cal = Cal.combine(CalEntry);

%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-12-01T00:00:00Z/2020-12-14T23:59:59Z'),...
  repmat([0.6834  3.8871],2,1));
Cal = Cal.combine(CalEntry);

%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-12-15T00:00:00Z/2020-12-31T23:59:59Z'),...
  repmat([0.5556  3.8249],2,1));
Cal = Cal.combine(CalEntry);

%-----------------------2021 ---------------------------------------------- 
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-01-01T00:00:00Z/2021-01-06T23:59:59Z'),...
  repmat([0.5905  4.0923],2,1)); 
Cal = Cal.combine(CalEntry);

%2
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-01-07T00:00:00Z/2021-01-12T05:49:59Z'),...
  repmat([0.6730  4.1837],2,1)); 
Cal = Cal.combine(CalEntry);


%3
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-01-12T05:50:00Z/2021-01-17T23:59:59Z'),...
  repmat([0.7462  4.5630],2,1)); 
Cal = Cal.combine(CalEntry);

%4
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-01-18T00:00:00Z/2021-02-16T23:59:59Z'),...
  repmat([0.3524  3.8309],2,1)); 
Cal = Cal.combine(CalEntry);


%5
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-02-17T00:00:00Z/2021-03-20T01:29:59Z'),...
  repmat([0.7651  3.4971],2,1)); 
Cal = Cal.combine(CalEntry);



%6
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-03-20T01:30:00Z/2021-03-22T19:29:59Z'),...
  repmat([0.6460 + 3.5047i   0.3899 + 3.7683i],2,1)); 


PSPintersection = 1.0297; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%   



%7
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-03-22T19:30:00Z/2021-04-04T03:59:59Z'),...
  repmat([0.7884 3.3714],2,1)); 
Cal = Cal.combine(CalEntry);

%8
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-04-04T04:00:00Z/2021-07-25T23:59:59Z'),...
  repmat([0.7125 + 3.0114i   0.4926 +  3.2371i],2,1)); 


PSPintersection = 1.050; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%  



%9
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-07-26T00:00:00Z/2021-08-02T23:59:59Z'),...
  repmat([0.7396 3.2209],2,1)); 
Cal = Cal.combine(CalEntry);

%10
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-08-03T00:00:00Z/2021-08-05T07:59:59Z'),...
  repmat([0.7694 3.3844],2,1)); 
Cal = Cal.combine(CalEntry);

%11
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-08-05T08:00:00Z/2021-08-08T23:59:59Z'),...
  repmat([0.6615 3.1782],2,1)); 
Cal = Cal.combine(CalEntry);

%12
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-08-09T00:00:00Z/2021-08-11T23:59:59Z'),...
  repmat([0.6802 3.2830],2,1)); 
Cal = Cal.combine(CalEntry);



%13
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-08-12T00:00:00Z/2021-08-29T19:59:59Z'),...
  repmat([0.7572 + 3.3979i  0.4088 + 3.6251i],2,1)); 
  

PSPintersection = 0.6522; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------% 



%14
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-08-29T18:00:00Z/2021-08-31T23:59:59Z'),...
  repmat([0.6898 + 3.7050i  0.4299 + 3.6844i],2,1)); 
  

PSPintersection = -0.0793; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------% 


%15
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-09-01T00:00:00Z/2021-09-04T00:59:59Z'),...
  repmat([0.6245 + 3.5472i  0.2589 + 3.8734i],2,1)); 
  

PSPintersection = 0.8925; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------% 

%16
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-09-04T01:00:00Z/2021-09-11T23:59:59Z'),...
  repmat([0.6916  3.5999],2,1)); 
Cal = Cal.combine(CalEntry);


%17
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-09-12T00:00:00Z/2021-09-27T22:05:45Z'),...
  repmat([0.5830 + 3.6131i  0.4904 + 3.5927i],2,1)); 


PSPintersection = 0.6390; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%

%18
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-09-27T22:05:46Z/2021-10-05T01:59:59Z'),...
  repmat([0.7342 + 3.5247i  0.4475 + 3.4873i],2,1)); 


PSPintersection = -0.1302; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%

%19
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-10-05T02:00:00Z/2021-10-09T23:59:59Z'),...
  repmat([0.7239  3.3958],2,1)); 
Cal = Cal.combine(CalEntry);


%20
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-10-10T00:00:00Z/2021-10-14T23:59:59Z'),...
  repmat([0.8061  3.2660],2,1)); 
Cal = Cal.combine(CalEntry);

%21
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-10-15T00:00:00Z/2021-10-15T23:59:59Z'),...
  repmat([0.6726  3.2377],2,1)); 
Cal = Cal.combine(CalEntry);

%22
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-10-16T00:00:00Z/2021-10-23T23:59:59Z'),...
  repmat([0.8444 + 3.1733i  0.5509 + 3.3010i],2,1)); 


PSPintersection = 0.550; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%

%23
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-10-24T00:00:00Z/2021-11-05T23:59:59Z'),...
  repmat([0.6579 + 3.1300i  0.4599 + 3.2233i],2,1)); 


PSPintersection = 0.4709; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%

%24
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-11-06T00:00:00Z/2021-11-20T23:59:59Z'),...
  repmat([0.7547  2.8010],2,1)); 
Cal = Cal.combine(CalEntry);


%25
%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-11-21T00:00:00Z/2021-12-26T23:59:59Z'),...
  repmat([0.6837 + 2.7283i  0.5384 + 2.8575i],2,1)); 


PSPintersection = 0.710; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%

%26
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-12-27T00:00:00Z/2021-12-31T23:59:59Z'),...
  repmat([0.5467  2.7777],2,1)); 
Cal = Cal.combine(CalEntry);


%% calibrate
CalR = Cal.resample(PSP);
NeScp = PSP; 


NeScp.data = exp(CalR.x.data.*NeScp.data + CalR.y.data);


timeOutsideInterval = irf_time('2021-12-31T23:59:59Z','utc>ttns');
NeScp.data(NeScp.time.epoch > timeOutsideInterval)= NaN;
    


NeScp.name = 'NeScp';
NeScp.units = 'cm^-3';
NeScp.siConversion = 'cm^-3>1e6*m^-3';
NeScp.userData = '';

end 


function [C] = TwofitCalibration(PSPint,y_eq,CalData)
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
        
        
        %Merge both fits and NaNs to keep the same lenght as the input
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


