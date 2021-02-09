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
codeVerStr = '2021-02-09T17:13:00';



% Based on data from 2020-04-07
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-03-08T00:00:00Z/2020-05-18T04:05:54Z'),...
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

% data until Sept 4
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-08-11T21:27:03Z/2020-09-04T23:59:59Z'),...
  repmat([0.2863  1.6844],2,1));

Cal = Cal.combine(CalEntry);

% data for Sept 5
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-09-05T00:00:00Z/2020-09-13T23:59:59Z'),...
  repmat([0.3104  1.5638],2,1));
Cal = Cal.combine(CalEntry);

%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-09-14T00:00:00Z/2020-11-02T23:59:59Z'),...
  repmat([0.2921  1.4582],2,1));
Cal = Cal.combine(CalEntry);

%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-11-03T00:00:00Z/2020-11-30T23:59:59Z'),...
  repmat([0.3199  1.5274],2,1));
Cal = Cal.combine(CalEntry);

%
CalEntry = irf.ts_vec_xy(...
  irf.tint('2020-12-01T00:00:00Z/2020-12-14T23:59:59Z'),...
  repmat([0.2956  1.6881],2,1));
Cal = Cal.combine(CalEntry);

%% calibrate
CalR = Cal.resample(PSP);
NeScp = PSP; 

NeScp.data = 10.^(CalR.x.data.*NeScp.data +CalR.y.data);

NeScp.name = 'NeScp';
NeScp.units = 'cm^-3';
NeScp.siConversion = 'cm^-3>1e6*m^-3';
NeScp.userData = '';
