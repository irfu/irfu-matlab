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
codeVerStr = '2022-04-08T16:45:00';


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


CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-02-02T00:00:00Z/2021-04-04T23:59:50Z'),...
  repmat([0.8000 + 3.4468i   0.3906 + 3.7661i],2,1)); 
Cal = Cal.combine(CalEntry);


CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-04-05T00:00:00Z/2021-07-27T23:59:59Z'),...
  repmat([0.7092 3.0440],2,1)); 
Cal = Cal.combine(CalEntry);


CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-07-28T00:00:00Z/2021-09-04T23:59:50Z'),...
  repmat([0.7812 + 3.3793i  0.3953 + 3.6551i],2,1)); 
Cal = Cal.combine(CalEntry);


CalEntry = irf.ts_vec_xy(...
  irf.tint('2021-09-05T00:00:00Z/2021-11-22T23:59:59Z'),...
  repmat([0.5882  3.4788],2,1)); 
Cal = Cal.combine(CalEntry);

%% calibrate
CalR = Cal.resample(PSP);
NeScp = PSP; 


NeScp.data = exp(real(CalR.x.data).*NeScp.data + real(CalR.y.data));

%Calibration intervals with 2 fits (i.e when CalR is complex)
fit2_time = Cal.time(imag(Cal.data(:,1)) ~= 0);

checkInterval_1 = PSP.tlim(fit2_time(1:2));

if ~isempty(checkInterval_1)
%=========================================================================
%2 fit- interval # 1 --> '2021-02-02T00:00:00Z/2021-04-04T23:59:59Z'


    %Intersection of the 2 fits
    y_eq = 0.7814;
    
    %Time interval for the fit with 2 coefficients #1
    Tstart = find(NeScp.time==checkInterval_1.time(1)); 
    Tend = find(NeScp.time==checkInterval_1.time(end));

    
    %Create a vector with time values with fit 1 (PSP1) and a vector with values
    %with fit 2 (PSP2) based on the intersection point y_eq
    PSP1.time = PSP.tlim(fit2_time(1:2)).time(PSP.tlim(fit2_time(1:2)).data<y_eq);
    PSP2.time = PSP.tlim(fit2_time(1:2)).time(PSP.tlim(fit2_time(1:2)).data>=y_eq);
   

    %Create a vector with data values with fit 1 (PSP1) and a vector with values
    %with fit 2 (PSP2) based on the intersection point y_eq
    PSP1.data = PSP.tlim(fit2_time(1:2)).data(PSP.tlim(fit2_time(1:2)).data<y_eq);
    PSP2.data = PSP.tlim(fit2_time(1:2)).data(PSP.tlim(fit2_time(1:2)).data>=y_eq);
    
    %Create a vector with NAN values and their times that will be incorporated in the final
    %TSeries to maintain the same length as the input
    PSPnan.time = PSP.tlim(fit2_time(1:2)).time(isnan(PSP.tlim(fit2_time(1:2)).data));
    PSPnan.data = PSP.tlim(fit2_time(1:2)).data(isnan(PSP.tlim(fit2_time(1:2)).data));
    
    
    %Create a vector with the coefficient values for the 2 fits based on
    %the intersection point y_eq
    CalR1.data = CalR.tlim(fit2_time(1:2)).x.data(PSP.tlim(fit2_time(1:2)).data<y_eq);
    CalR2.data = CalR.tlim(fit2_time(1:2)).y.data(PSP.tlim(fit2_time(1:2)).data>=y_eq);

    
    %Convert the vectors to TSeries
    NeScp_1 = TSeries(PSP1.time,PSP1.data);
    NeScp_2 = TSeries(PSP2.time,PSP2.data);
    NeScp_nan = TSeries(PSPnan.time,PSPnan.data);
    CalR1 = TSeries(PSP1.time,CalR1.data);
    CalR2 = TSeries(PSP2.time,CalR2.data);


    %Perform the 2 fits (above and below the intersection point)
    NeScp_1.data = exp(real(CalR1.data).*NeScp_1.data +imag(CalR1.data));
    NeScp_2.data = exp(real(CalR2.data).*NeScp_2.data +imag(CalR2.data));

    
    %Check if the vectors are empty (i.e. 2 fits are not needed)
    if ~isempty(NeScp_1) && ~isempty(NeScp_2)
        NeScp_2fits = NeScp_2.combine(NeScp_1);
    elseif ~isempty(NeScp_1) && isempty(NeScp_2)
        NeScp_2fits = NeScp_1;
    elseif  ~isempty(NeScp_2) && isempty(NeScp_1)
       NeScp_2fits = NeScp_2;
    elseif isempty(NeScp_1) && isempty(NeScp_2)
        NeScp_2fits = NeScp(Tstart:Tend);
    end
    
    %Check NANs to keep the same length of NeScp as the original
    if ~isempty(NeScp_nan)
        NeScp_2fits = NeScp_2fits.combine(NeScp_nan);
    end
    
    %Incorporate the fitted values with 2 coefficients to the rest of the
    %data
    NeScp.data(Tstart:Tend) = NeScp_2fits.data;
    
end


checkInterval_2 = PSP.tlim(fit2_time(3:4));

if ~isempty(checkInterval_2)
%==============================================================================
    clear PSP1 PSP2 PSPnan CalR1 CalR2 NeScp_1 NeScp_2 y_eq NeScp_2fits Tstart Tend
%==============================================================================


%2 fit- interval # 2 --> '2021-07-28T00:00:00Z/2021-09-04T23:59:59Z'

    %Intersection of the 2 fits
    y_eq = 1.7180;
    
    %Time interval for the fit with 2 coefficients #2
    Tstart = find(NeScp.time==checkInterval_2.time(1)); 
    Tend = find(NeScp.time==checkInterval_2.time(end));

    
    %Create a vector with time values with fit 1 (PSP1) and a vector with values
    %with fit 2 (PSP2) based on the intersection point y_eq
    PSP1.time = PSP.tlim(fit2_time(3:4)).time(PSP.tlim(fit2_time(3:4)).data<y_eq);
    PSP2.time = PSP.tlim(fit2_time(3:4)).time(PSP.tlim(fit2_time(3:4)).data>=y_eq);
    
    %Create a vector with data values with fit 1 (PSP1) and a vector with values
    %with fit 2 (PSP2) based on the intersection point y_eq
    PSP1.data = PSP.tlim(fit2_time(3:4)).data(PSP.tlim(fit2_time(3:4)).data<y_eq);
    PSP2.data = PSP.tlim(fit2_time(3:4)).data(PSP.tlim(fit2_time(3:4)).data>=y_eq);
    
    %Create a vector with NAN values and their times that will be incorporated in the final
    %TSeries to maintain the same length as the input
    PSPnan.data = PSP.tlim(fit2_time(3:4)).data(isnan(PSP.tlim(fit2_time(3:4)).data));
    PSPnan.time = PSP.tlim(fit2_time(3:4)).time(isnan(PSP.tlim(fit2_time(3:4)).data));

    %Create a vector with the coefficient values for the 2 fits based on
    %the intersection point y_eq
    CalR1.data = CalR.tlim(fit2_time(3:4)).x.data(PSP.tlim(fit2_time(3:4)).data<y_eq);
    CalR2.data = CalR.tlim(fit2_time(3:4)).y.data(PSP.tlim(fit2_time(3:4)).data>=y_eq);

    %Convert the vectors to TSeries
    NeScp_1 = TSeries(PSP1.time,PSP1.data);
    NeScp_2 = TSeries(PSP2.time,PSP2.data);
    NeScp_nan = TSeries(PSPnan.time,PSPnan.data);
    CalR1 = TSeries(PSP1.time,CalR1.data);
    CalR2 = TSeries(PSP2.time,CalR2.data);


    %Perform the 2 fits (above and below the intersection point)
    NeScp_1.data = exp(real(CalR1.data).*NeScp_1.data +imag(CalR1.data));
    NeScp_2.data = exp(real(CalR2.data).*NeScp_2.data +imag(CalR2.data));


    %Check if the vectors are empty (i.e. 2 fits are not needed)
    if ~isempty(NeScp_1) && ~isempty(NeScp_2)
        NeScp_2fits = NeScp_2.combine(NeScp_1);
    elseif ~isempty(NeScp_1) && isempty(NeScp_2)
        NeScp_2fits = NeScp_1;
    elseif  ~isempty(NeScp_2) && isempty(NeScp_1)
       NeScp_2fits = NeScp_2;
    elseif isempty(NeScp_1) && isempty(NeScp_2)
        NeScp_2fits = NeScp(Tstart:Tend);
    end
    
    %Check NANs to keep the same length of NeScp as the original
    if ~isempty(NeScp_nan)
        NeScp_2fits = NeScp_2fits.combine(NeScp_nan);
    end
    
    %Incorporate the fitted values with 2 coefficients to the rest of the
    %data
    NeScp.data(Tstart:Tend) = NeScp_2fits.data;
%================================================================
end

timeOutsideInterval = irf_time('2021-11-22T23:59:59Z','utc>ttns');
NeScp.data(NeScp.time.epoch > timeOutsideInterval)= NaN;
    


NeScp.name = 'NeScp';
NeScp.units = 'cm^-3';
NeScp.siConversion = 'cm^-3>1e6*m^-3';
NeScp.userData = '';
