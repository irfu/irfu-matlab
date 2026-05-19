function [DCE_SRF_out, PSP_out, ScPot_out, codeVerStr, matVerStr] = vdccal(VDC_inp, EDC_inp, calFilename)
%SOLO.VDCCAL  Calibrate VDC to get DC E and PSP
%
%    [DCE_SRF,PSP,ScPot,codeVerStr,matVerStr] = solo.vdccal(VDC,EDC,calfilename)
%
%
% ARGUMENTS
% =========
% VDC_inp, EDC_inp
%       TSeries objects with data from L2 CWF file(s).
%       NOTE: BICAS calls this function but sets samples to NaN for timestamps
%             for which a "synthetic" version of L2 zVariable QUALITY_FLAG is
%             (strictly) lower than the value of BICAS setting
%             PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN. This value is set to "2"
%             by default.
%       /Erik P G Johansson, 2026-05-18
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
%       Must be on the form of a human-readable UTC timestamp string which
%       conforms to regular expression bicas.proc.L2L3.ext.CODE_VER_STR_REGEXP.
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



% Assert that all calibration data arrays use the same timestamps.
assert(all((a.d23.time == a.K123.time) & (a.d23.time == a.k23.time)))

% Extract timestamps for beginning and end of calibration data, so that one can
% implement special handling for time for which there is no calibration data.
timeCalibrationDataBegin = a.d23.time(1);
timeCalibrationDataEnd   = a.d23.time(end);



%=============================================================
% Set return values that represent the version of calibration
%=============================================================
% Version of the function (not .mat file)
% ---------------------------------------
% This string is used by BICAS to set a CDF global attribute in official
% datasets for traceability.
% NOTE: This value is meant to be be UPDATED BY HAND, not by an automatic
%       timestamp, so that a constant value represents the same
%       function/algorithm.
codeVerStr = '2026-05-19T14:15:00';
% Version of the .mat file
% ------------------------
% This string is used by BICAS to set a CDF global attribute in official
% datasets for traceability.
% Using filename, or at least for now.
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
calTint      = irf.tint(a.d23.time(1), a.d23.time(end));
% We only care about discontinuities inside the calibration interval.
discontTimes = EpochTT( ...
  allDiscontTimes.epoch(allDiscontTimes.epoch <= calTint.epoch(2) )...
  );

mainTint      = irf.tint(VDC_inp.time(1), VDC_inp.time(end));
sub_int_times = EpochTT(solo.split_tint(mainTint, discontTimes));

% Predefine empty output variables:
DCE_SRF_out = irf.ts_vec_xyz(EpochTT([]), double.empty(0, 3));
PSP_out     = irf.ts_scalar( EpochTT([]), double.empty(0, 1));
ScPot_out   = irf.ts_scalar( EpochTT([]), double.empty(0, 1));



% Perform calibration on each subinterval separately (if any probe-to-spacecraft
% potential discontinuities are present).

for iSub = 1:length(sub_int_times)-1
  % =========================================================
  % Identify time interval, defined by listed discontinuities
  % =========================================================

  subTint = sub_int_times(iSub:iSub+1);
  % Find the closest discontinuities.
  prev_discont = EpochTT(max(discontTimes.epoch(subTint(1).epoch   >= discontTimes.epoch)));
  next_discont = EpochTT(min(discontTimes.epoch(subTint(end).epoch <= discontTimes.epoch)));

  % -----------------------------------
  % Optionally modify the time interval
  % -----------------------------------
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



  % =====================================================================
  % Prepare for calibration: Same timestamps for VDC and calibration data
  % =====================================================================

  VDC = VDC_inp.tlim(subTint);

  % VDC indices for which there is no calibration data.
  bNoCalibrationData = ...
    (VDC.time < timeCalibrationDataBegin) | ...
    (VDC.time > timeCalibrationDataEnd);

  % -------------------------------------------------
  % Resample calibration parameters to VDC timestamps
  % -------------------------------------------------
  % NOTE: Empirically, .resample() uses linear interpolation.
  % IMPORTANT NOTE: .resample() somehow extrapolates data to outside of the time
  % interval for which there is data. ==> Generates new "calibration data"
  % values (not NaN) also outside of time interval for which there is actual
  % valid calibration data!
  d23R  = a.d23.tlim( subTint).resample(VDC);
  k23R  = a.k23.tlim( subTint).resample(VDC);
  K123R = a.K123.tlim(subTint).resample(VDC);
  % Set (resampled) calibration data to NaN for timestamps outside of the time
  % interval for which there is actual calibration data. One must do this to
  % prevent wildly extrapolating calibration data to far outside the time
  % interval covered by actual calibration data, which leads to calibrating data
  % which can not (yet) be calibrated.
  d23R.data( bNoCalibrationData) = NaN;
  k23R.data( bNoCalibrationData) = NaN;
  K123R.data(bNoCalibrationData) = NaN;



  % =====================================
  % Calibrate: Derive PSP, ScPot, DCE_SRF
  % =====================================

  V1         = double(VDC.x.data);
  % Remove potential offset between probes 2 & 3.
  V2_scaled  = double(VDC.y.data) .* k23R.data + double(d23R.data);
  V3         = double(VDC.z.data);
  % V23: Corresponds to a measurement point between antennas 2 & 3.
  V23        = (V2_scaled+V3) / 2;

  V23_scaled = V23 .* K123R.data(:,1) + K123R.data(:,2); % Correcting V23 to V1.

  % ----------
  % Derive PSP
  % ----------
  % Assume all probe data available: Compute PSP from corrected quantities.
  PSP = irf.ts_scalar(VDC.time, (V23_scaled + V1)/2);
  % Use alternate "calculation" using only ANT1 for some timestamps.
  % NOTE: This calculation is independent of calibration data. Should therefore
  %       not use NaN when bNoCalibrationData==true.
  bPspOnlyUsesAnt1 = isnan(VDC.y.data) | isnan(VDC.z.data);
  bPspOnlyUsesAnt1 = bPspOnlyUsesAnt1  | bNoCalibrationData;
  PSP.data(bPspOnlyUsesAnt1) = VDC.x.data(bPspOnlyUsesAnt1);
  %
  PSP.units = 'V';
  PSP_out   = PSP_out.combine(PSP);

  % -----------------------------
  % Derive ScPot: Function of PSP
  % -----------------------------
  % XXX: these are just ad hoc numbers.
  PLASMA_POT   = 1.5;
  SHORT_FACTOR = 2.5;
  %
  ScPot = irf.ts_scalar(VDC.time, -(PSP.data-PLASMA_POT)*SHORT_FACTOR);
  ScPot.units = PSP.units;
  ScPot_out   = ScPot_out.combine(ScPot);

  % --------------
  % Derive DCE_SRF
  % --------------
  % Ey_SRF = V3 - V2, 6.99 - 1/2 of distance between the antennas
  V_d23  = V2_scaled-V3;    % Fixed V2-V3.
  Ey_SRF = -V_d23*1e3/6.99;

  % Here we use the effective antenna length of 11.2 m, which correponds to
  % the distance between the center of ANT1 and a symmetric antenna on the
  % other side having voltage V23_scaled (see above).

  % For non-scaled V23, the effective length would be L_123 = 6.97m, as shown
  % in Steinvall et al., 2021.
  Ez_SRF = (V23_scaled - V1)*1e3/11.2;

  % NOTE: Setting X component to zero (except when Ey_SRF is NaN).
  % NOTE: Ey_SRF may contain NaN. Therefore Ey_SRF*0 != zeros(size(Ey_SRF)).
  %       (Bug?!)
  DCE_SRF = irf.ts_vec_xyz(VDC.time, [Ey_SRF*0, Ey_SRF, Ez_SRF]);
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

end % function
