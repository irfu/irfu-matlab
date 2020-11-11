function [DCE_SRF,PSP,ScPot] = vdccal(VDC,EDC)
%SOLO.VDCCAL  Calibrate VDC to get DC E and PSP
%
%    [DCE_SRF,PSP,ScPot] = solo.vdccal(VDC,EDC)
%
% Inputs: VDC,EDC from L2 CWF files
%
% Outputs:
%   DCE_SRF - DC electric field in SRF (Ex=0)
%   PSP     - probe-to-spacecraft potential
%   ScPot   - spacecraft potential (PRELIMINARY PROXY)
%
% Loads d23K123.mat file produced by solo.correlate_probes_batch (script)
% NOTE: d23K123.mat needs to be updated before processing the new data.
%
% NOTE: This function is used by BICAS for producing official datasets.

a = load('d23K123.mat');

d23R = a.d23.resample(VDC);
K123R = a.K123.resample(VDC);

V2corr = double(VDC.y.data) -double(d23R.data);

V23 = (V2corr + double(VDC.z.data))/2; % (V2 + V3) /2
V23corr = (V23.*K123R.data(:,1) + K123R.data(:,2));

PSP = irf.ts_scalar(VDC.time,(V23corr + double(VDC.x.data))/2);
PSP.units = 'V';

PLASMA_POT = 1.5; SHORT_FACTOR = 2.5; % XXX: these are just adhoc numbers

ScPot = irf.ts_scalar(VDC.time, -(PSP.data-PLASMA_POT)*SHORT_FACTOR);
ScPot.units = PSP.units;

% Ey_SRF = V3 - V2, 6.99 - 1/2 of distance between the antennas
E23 = double(EDC.z.data) -double(d23R.data);
Ey_SRF = -E23*1e3/6.99;

% Ez_SRF = V23 - V1
%E12 = double(EDC.x.data) - PSP.data*0.1269; % correct for common mode
%Ez_SRF = -(E12 + E23/2) *1e3/6.97;
Ez_SRF = V23corr - double(VDC.x.data);	% V23 - V1	
Ez_SRF = Ez_SRF*1e3/6.97;

DCE_SRF = irf.ts_vec_xyz(VDC.time,[Ey_SRF*0 Ey_SRF Ez_SRF]);
DCE_SRF.units = 'mV/m';
DCE_SRF.coordinateSystem = 'SRF';

