function [DCE_SRF,PSP,SCPOT] = vdccal(VDC)
%VDCCAL  Calibrate VDC to get DC E and PSP
%
%    [DCE_SRF,PSP,SCPOT] = vdccal(VDC)
%
% NOTE: This function is used by BICAS for producing official datasets.
%
% d23K123.mat file produced by solo.correlate_probes_batch (script)
% d23K123.mat needs to be updated before processing the new data.

a = load('d23K123.mat');

d23R = a.d23.resample(VDC);
K123R = a.K123.resample(VDC);

V2corr = double(VDC.y.data) -double(d23R.data);

V23 = (V2corr + double(VDC.z.data))/2; % (V2 + V3) /2
V23corr = (V23.*K123R.data(:,1) + K123R.data(:,2));

E23 = double(VDC.z.data) - V2corr;
E23 = E23*1e3/6;

E123 = V23corr - double(VDC.x.data);
E123 = E123*1e3/6;

DCE_SRF = irf.ts_vec_xyz(VDC.time,[E23*0 -E23 E123]);
DCE_SRF.units = 'mV/m';
DCE_SRF.coordinateSystem = 'SRF';

PSP = irf.ts_scalar(VDC.time,(V23corr + double(VDC.x.data))/2);
PSP.units = 'V';

SCPOT = irf.ts_vec_xyz(VDC.time, PSP.data*[-1,-1,-1] + [0,0,0]);
SCPOT.units = PSP.units;