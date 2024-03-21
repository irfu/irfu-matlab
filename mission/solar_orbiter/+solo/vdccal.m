function [DCE_SRF_out, PSP_out, ScPot_out, codeVerStr,matVerStr] = vdccal(VDC_inp, EDC_inp, calFilename)
%SOLO.VDCCAL  Calibrate VDC to get DC E and PSP
%
%    [DCE_SRF,PSP,ScPot,codeVerStr,matVerStr] = solo.vdccal(VDC,EDC,calfilename)
%
%
% ARGUMENTS
% =========
% VDC_inp, EDC_inp
%       TSeries objects with data from L2 CWF file(s).
%       NOTE: 2023-10-04: BICAS calls this function but sets samples to NaN for
%       timestamps for which QUALITY_FLAG are (strictly) lower than BICAS
%       setting PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN which is set to "2" by
%       default.
%       VDC_inp: x, y, z = V1,  V2,  V3
%       EDC_inp: x, y, z = V12, V13, V23
% calFilename
%       The name of the calibration file (.mat) one wishes to
%       use. If empty, then code will use BICAS's official calibration file
%       produced by solo.correlate_probes_batch (script) (?) and that is used
%       for producing official datasets.
%
%
% RETURN VALUES
% =============
% DCE_SRF
%       DC electric field in SRF (Ex=0).
% PSP
%       Probe-to-spacecraft potential
% ScPot
%       Spacecraft potential (PRELIMINARY PROXY)
% codeVerStr
%       Date-formatted version string for the code ("algorithm") that implements
%       the function itself. This is used by BICAS for setting the relevant CDF
%       global attribute to indicate the version of the algorithm used to
%       produce a particular dataset (CDF file).
%       Note: The version string does not cover the calibration .mat file.
% matVerStr
%       Version string representing the calibration .mat file. Currently
%       consists of the calibration file filename which in turn contains the
%       relevant versioning date. This is used by BICAS for setting the relevant
%       CDF global attribute to indicate the calibration data file used to
%       produce a particular dataset (CDF file).
%
%
% NOTE: Calibration .mat file only has a certain time coverage and therefore
%       needs to be updated before processing new data outside the already
%       covered time interval.
% NOTE: This function is used by BICAS for producing official L3 datasets. It
%       must therefore have an interface that is compatible with BICAS.


% Timestamp after which all data should be regarded as having only one probe of
% data (VDC1) available.
% -----------------------------------------------------------------------------
% NOTE: This is to mitigate against a long period of bad data with a lot of
% saturated DC diffs, beginning somewhere around Dec 2022 and that has not yet
% ended as of 2023-08-17. There will hopefully eventually be an associated end
% date to this time period of bad data. This special treatment may or may not be
% a temporary measure.
% PROPOSAL: Have BICAS exclude saturated VDC diffs being sent to this function
% in the first place. BICAS could exclude VDC diffs based on
% (not-yet-implemented) quality bits or NSO table.
TIME_PSP_BEGIN_SINGLE_PROBE = EpochTT('2022-12-15T00:00:00.000000000Z');



% Normalize "calFilename": Always contain filename.
if isempty(calFilename)
  % Caller did not specify calibration file.
  % IMPORTANT: USES CALIBRATION FILE THAT IS USED BY BICAS FOR PRODUCING
  % OFFICIAL DATASETS.
  calFilename = 'd23K123_20230707.mat'; % parameters up to end of 2023-05-27
else
  % Caller specified calibration file. Useful for debugging/testing new
  % calibrations.

  % (Do nothing.)
end
a = load(calFilename);



%=============================================================
% Set return values that represent the version of calibration
%=============================================================
% Version of the function (not .mat file).
% NOTE: This value is meant to be be UPDATED BY HAND, not by an automatic
% timestamp, so that a constant value represents the same function/algorithm.
codeVerStr = '2022-12-06T13:23:14';
% Version of the .mat file. Using filename, or at least for now.
% This string is used by BICAS to set a CDF global attribute in official
% datasets for traceability.
[~, basename, suffix] = fileparts(calFilename);
matVerStr  = [basename, suffix];   % Only use filename, not entire path.

%=============================================================================
% Find data points/CDF records for which only V1_DC is available
%
% NOTE: solo.vdccal() only uses DC data (not AC). For mux=0, there are two
% cases:
% (1) all single probes (DC) available
% (2) only probe 1 (DC) available.
% NOTE: Only works for mux=0,2,3,4 (not mux=1).
% NOTE: Currently ignores EDC_inp argument.
%==============================================================================

% Load Probe-potential discontinuities.
allDiscontTimes = solo.ProbePotDiscontinuities;
% Time interval defined by the calibration file.
calTint      = irf.tint(a.d23.time(1),a.d23.time(end));
% We only care about discontinuities inside the calibration interval.
discontTimes = EpochTT(allDiscontTimes.epoch(allDiscontTimes.epoch<=calTint.epoch(2)));

mainTint      = irf.tint(VDC_inp.time(1),VDC_inp.time(end));
sub_int_times = EpochTT(solo.split_tint(mainTint,discontTimes));

% Predefine output variables:
DCE_SRF_out = irf.ts_vec_xyz(EpochTT([]),double.empty(0,3));
PSP_out     = irf.ts_scalar(EpochTT([]),[]);
ScPot_out   = irf.ts_scalar(EpochTT([]),[]);

% Perform calibration on each subinterval separately (if any probe-to-spacecraft
% potential discontinuities are present).

for iSub = 1:length(sub_int_times)-1

  subTint = sub_int_times(iSub:iSub+1);
  % Find the closest discontinuities.
  prev_discont = EpochTT(max(discontTimes.epoch(subTint(1).epoch   >= discontTimes.epoch)));
  next_discont = EpochTT(min(discontTimes.epoch(subTint(end).epoch <= discontTimes.epoch)));

  % ======================================
  % Optionally modify the time subinterval
  % ======================================
  % Extend the time interval to the closest discontinuities (this is to avoid
  % problems that occur if VDC contains very little data).
  if ~isempty(prev_discont) && ~isempty(next_discont)
    subTint = irf.tint(prev_discont, next_discont);
  elseif isempty(prev_discont)
    % If there are no discontinuities BEFORE the specified time,
    % increase interval by 2 days before.
    subTint = irf.tint(subTint(1)+(-2*24*60*60), next_discont);
  elseif isempty(next_discont)
    % If there are no discontinuities AFTER the specified time,
    % increase interval by 2 days after.
    subTint = irf.tint(prev_discont, subTint(end)+(2*24*60*60));
  end

  %%
  % =======================
  % Prepare for calibration
  % =======================
  VDC = VDC_inp.tlim(subTint);

  % Indices/samples for which data should be treated as single probe.
  bSingleProbe = isnan(VDC.y.data) & isnan(VDC.z.data);
  bSingleProbe = bSingleProbe | (VDC.time > TIME_PSP_BEGIN_SINGLE_PROBE);

  % Resample calibration parameters
  d23R  = a.d23.tlim(subTint).resample(VDC);
  k23R  = a.k23.tlim(subTint).resample(VDC);
  K123R = a.K123.tlim(subTint).resample(VDC);

  % =================
  % Begin calibration
  % =================

  V1         = double(VDC.x.data);
  % Remove potential offset between probes 2 & 3.
  V2_scaled  = double(VDC.y.data).*k23R.data + double(d23R.data);
  V3         = double(VDC.z.data);
  % V23: Corresponds to a measurement point between antennas 2 & 3.
  V23        = (V2_scaled+V3)/2;

  V23_scaled = (V23.*K123R.data(:,1) + K123R.data(:,2)); % Correcting V23 to V1.

  % Assume all probe data available: Compute PSP from corrected quantities.
  PSP = irf.ts_scalar(VDC.time, (V23_scaled + V1)/2);
  % Single-probe data: Use alternate, simpler "calculation" for some
  %                    timestamps.
  PSP.data(bSingleProbe) = VDC.x.data(bSingleProbe);
  PSP.units = 'V';
  PSP_out   = PSP_out.combine(PSP);

  % XXX: these are just ad hoc numbers.
  PLASMA_POT   = 1.5;
  SHORT_FACTOR = 2.5;

  ScPot = irf.ts_scalar(VDC.time, -(PSP.data-PLASMA_POT)*SHORT_FACTOR);
  ScPot.units = PSP.units;
  ScPot_out   = ScPot_out.combine(ScPot);

  % Ey_SRF = V3 - V2, 6.99 - 1/2 of distance between the antennas
  V_d23  = V2_scaled-V3;    % Fixed V2-V3.
  Ey_SRF = -V_d23*1e3/6.99;

  % Here we use the effective antenna length of 11.2 m, which correponds to
  % the distance between the center of ANT1 and a symmetric antenna on the
  % other side having voltage V23_scaled (see above).

  % For non-scaled V23, the effective length would be L_123 = 6.97m, as shown
  % in Steinvall et al., 2021.
  Ez_SRF = (V23_scaled - V1)*1e3/11.2;

  % NOTE: Ey_SRF may contain NaN. Therefore Ey_SRF*0 != zeros(size(Ey_SRF)).
  % (Bug?!)
  DCE_SRF = irf.ts_vec_xyz(VDC.time, [Ey_SRF*0 Ey_SRF Ez_SRF]);
  DCE_SRF.units            = 'mV/m';
  DCE_SRF.coordinateSystem = 'SRF';

  DCE_SRF_out = DCE_SRF_out.combine(DCE_SRF);

end % for

% Specify units and coordinate system for the variables that are actually
% returned from the function.
% -----------------------------------------------------------------------
% NOTE: Set .units and .coordinateSystem explicitly since
% (1) TSeries.combine(), which fill the objects with science data, will not
%     copy those values, and
% (2) they can not be copied from PSP, ScPot, DCE_SRF since they will not be set
%     if there are zero sub-time intervals (zero loop iterations).
DCE_SRF_out.units            = 'mV/m';
DCE_SRF_out.coordinateSystem = 'SRF';
PSP_out.units                = 'V';
ScPot_out.units              = 'V';

end %function
