function NeScp = psp2ne(PSP)
%SOSO.PSP2NE  Convert probe-to-spacecraft potential to electron density
%
% NeScp = solo.psp2ne(PSP)
%
% Convert probe-to-spacecraft (PSP) potential to electron density (NeScp)
%
% The calibration is based on the RPW/QTN/Fpe data

% Calibraiton using plasma line 
% see Dropbox/Solar_Orbiter/Science data/InFlight Cal/Ncalpsp2ne_calibrate.m
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-04-07T00:00:00Z/2020-05-18T04:05:54Z'),...
  repmat([0.3835   1.4908],2,1));

Cal = CalEntry;

CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-05-18T04:05:55Z/2020-05-29T23:59:59Z'),...
  repmat([0.3539   1.9785],2,1));

Cal = Cal.combine(CalEntry);

% cal based up to 2020-07-05T23:59:59Z
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-05-30T00:00:00Z/2020-08-11T21:27:02Z'),...
  repmat([0.2260   1.9106],2,1));

Cal = Cal.combine(CalEntry);

CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-08-11T21:27:03Z/2020-08-30T23:59:59Z'),...
  repmat([0.3116  1.6966],2,1));

Cal = Cal.combine(CalEntry);

%% calibrate
CalR = Cal.resample(PSP);
NeScp = PSP; 

NeScp.data = 10.^(CalR.x.data.*NeScp.data +CalR.y.data);
