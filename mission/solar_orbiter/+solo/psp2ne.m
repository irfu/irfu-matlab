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
codeVerStr = '2022-04-12T15:31:00';


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

%2021 -- 
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-01-01T00:00:00Z/2021-02-01T23:59:59Z'),...
  repmat([0.8000  4.5233],2,1)); 
Cal = Cal.combine(CalEntry);




%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-02-02T00:00:00Z/2021-04-04T23:59:50Z'),...
  repmat([0.8000 + 3.4468i   0.3906 + 3.7661i],2,1)); 


PSPintersection = 0.7814; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------%   



%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-04-05T00:00:00Z/2021-07-27T23:59:59Z'),...
  repmat([0.7092 3.0440],2,1)); 
Cal = Cal.combine(CalEntry);





%======================2 fits=========================%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-07-28T00:00:00Z/2021-09-04T23:59:50Z'),...
  repmat([0.7812 + 3.3793i  0.3953 + 3.6551i],2,1)); 
  

PSPintersection = 1.7180; %Intersection between 2 fits
checkInterval = PSP.tlim(CalEntry.time); %PSP data inside cal. interval

if ~isempty(checkInterval)
    [CalEntry] = TwofitCalibration(checkInterval,PSPintersection,CalEntry);
else
    CalEntry.data(1:end,2) = imag(CalEntry.x.data);
    CalEntry.data(1:end,1) = real(CalEntry.x.data);
end

Cal = Cal.combine(CalEntry);
%--------------------------------------------------------% 




%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-09-05T00:00:00Z/2021-11-22T23:59:59Z'),...
  repmat([0.5882  3.4788],2,1)); 
Cal = Cal.combine(CalEntry);


%% calibrate
CalR = Cal.resample(PSP);
NeScp = PSP; 


NeScp.data = exp(CalR.x.data.*NeScp.data + CalR.y.data);


timeOutsideInterval = irf_time('2021-11-22T23:59:59Z','utc>ttns');
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
        %If statements added for robusteness. An error will be send if 
        %combining empty objects.
        if ~isempty(C1) && ~isempty(C2)
            C = C1.combine(C2);
        elseif ~isempty(C1) && isempty(C2)
            C = C1;
        elseif ~isempty(C2) && isempty(C1)
           C = C2;
        elseif isempty(C1) && isempty(C2) && isempty(Cnan)
           C = CalData;
           irf.log('critical', 'no data at all ?!?')
        end
        
        if ~isempty(Cnan)
            C = C.combine(Cnan);
        end
    
    
end


