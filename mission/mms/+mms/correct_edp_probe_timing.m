function Vcorrected = correct_edp_probe_timing(Vorg)
% CORRECT_EDP_PROBE_TIMING corrects for the channel delays not accounted for in
% the MMS EDP processing. As described in the MMS EDP data products guide.
%
%  Vcorrected = CORRECT_EDP_PROBE_TIMING(Vorg);
%  Input:
%    Vorg       - TSeries created from L2 SCpot files, from the variable
%                 "mms#_edp_dcv_brst_l2" containing individual probe
%                 potentials.
%  Output:
%    Vcorrected - TSeries where the channel delay for each probe have been 
%                 accounted and corrected for.
%
% Note this function is only useful for Burst mode data. For the other
% telemetry modes (i.e. slow and fast) the channel delays are completly
% negligable and the interpolation and resampling applied here will have
% no effect other than possibly introduce numerical noise.

% Verify input
narginchk(1,1);
if(~isa(Vorg, 'TSeries'))
  errStr = ['Incorrect input class. Expected TSeries, got :', class(Vorg)];
  irf.log('critical', errStr); error(errStr);
end

% Show some comparision plots of results before and after (if true) or not
DEBUG = false;

% Reconstruct E12, E34, E56 as computed in MMS processing
timeOrig = Vorg.time;
E12 = TSeries(timeOrig, ( Vorg.data(:,1) - Vorg.data(:,2) )/0.120); 
E34 = TSeries(timeOrig, ( Vorg.data(:,3) - Vorg.data(:,4) )/0.120); 
E56 = TSeries(timeOrig, ( Vorg.data(:,5) - Vorg.data(:,6) )/0.0292);
% Correct the time tags to create individual time series
V1 = TSeries(timeOrig, Vorg.data(:,1));
V3 = TSeries(timeOrig + 7.629e-6, Vorg.data(:,3));
V5 = TSeries(timeOrig + 15.259e-6, Vorg.data(:,5));
E12.time = timeOrig + 26.703e-6;
E34.time = timeOrig + 30.518e-6;
E56.time = timeOrig + 34.332e-6;
% Resample all data to time tags of V1 (i.e. timeOrig).
V3 = V3.resample(timeOrig);
V5 = V5.resample(timeOrig);
E12 = E12.resample(timeOrig);
E34 = E34.resample(timeOrig);
E56 = E56.resample(timeOrig);
% Recompute individual probe potentials V2, V4, V6
V2 = V1 - E12 * 0.120;
V4 = V3 - E34 * 0.120;
V6 = V5 - E56 * 0.0292;
% Create the new TSeries with the corrected values
Vcorrected = irf.ts_scalar(timeOrig, [V1.data, V2.data, V3.data, V4.data, V5.data, V6.data]);

if(DEBUG)
  % Plot new values
  figure; %#ok<UNRCH>
  irf_plot({V1, V2, V3, V4, V5, V6});
  % Plot both old and new values
  figure;
  h = irf_plot({ ...
    TSeries(V1.time, [V1.data, Vorg.data(:,1)]), ...
    TSeries(V2.time, [V2.data, Vorg.data(:,2)]), ...
    TSeries(V3.time, [V3.data, Vorg.data(:,3)]), ...
    TSeries(V4.time, [V4.data, Vorg.data(:,4)]), ...
    TSeries(V5.time, [V5.data, Vorg.data(:,5)]), ...
    TSeries(V6.time, [V6.data, Vorg.data(:,6)])}, '.-');
  title(h(1), 'Corrected vs uncorrected channel delays.');
  legend(h(1), 'V1new', 'V1old');
  legend(h(2), 'V2new', 'V2old');
  legend(h(3), 'V3new', 'V3old');
  legend(h(4), 'V4new', 'V4old');
  legend(h(5), 'V5new', 'V5old');
  legend(h(6), 'V6new', 'V6old');
end

end
