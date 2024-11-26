% Try to autodetect sweeps.
%
% NOTE: This function is intended to be temporary functionality while waiting
% for the long-term solution which is to use the relevant bit in L1/L1R
% QUALITY_BITMASK for detecting sweeps. Since it is a (hopefully) short-term
% solution, the functionality is also not as sophisticated (and presumably
% accurate) or configurable as it could be.
%
% NOTE: There is some unimportant "imprecision" in the window algorithm
% which can be seen as an unimportant bug. Records which are BDM<>4 (i.e.
% definately not sweep) within a window with currents exceeding the
% thresholds (for BDM=4) are labelled as sweeps.
%
% NOTE: 2024-11-26: It *appear* as if L1 QUALITY_BITMASK in sample datasets
% (not production) contains bits for sweeps but without time margins before and
% after. Sweep bits in QUALITY_BITMASK has not been confirmed yet.
%
%
% ALGORITHM
% =========
% PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC specifies a timestamp.
% Records before timestamp:
%   BDN=4 <=> sweep
% Records after timestamp:
%   Among records with BDN=4, check sliding windows for the min-max
%   difference for BIAS HK's measured bias current (within the entire
%   window). For windows for which the min-max difference exceeds a
%   threshold, label entire window (where BDM=4) as sweeping.
%
%
% ARGUMENTS
% =========
% Bso
%       NOTE: PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS: If
%       greater than the number of CDF records/rows of data, then no
%       record will be labelled as sweeping.
%
function isSweepingFpa = autodetect_sweeps(hkTt2000, hkBdmFpa, hkBiasCurrentFpa, Bso)
% TODO-DEC: Does having argument and return value FPAs make sense?
%           Should caller convert?
% PROPOSAL: Rename (and negate) BDM_SWEEP_POSSIBLE --> BDM_SWEEP_IMPOSSIBLE
%   CON: Unintuitive for time period when sweep can be deduced frmo BDM.
% PROPOSAL: Rename "PTS" (unit) for PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS.
%   PRO: Unclear
%   PROPOSAL: "HK_CDF_RECORDS", HkCdfRecords
%     CON: Long
%
% PROPOSAL: Separate window margins for before and after window.
%   PRO: Margin after needs to be longer.
%   NOTE: Looking at sweep for 2024-06-21, margins should be maybe:
%       before sweep proper:  ~1-2 min
%       after sweep proper:   ~6-7 min
% PROPOSAL: Sweep detection algorithm which uses (and labels) the data gaps
%           before & after the sweep.
% PROPOSAL: Length of margins shouls be set in time, not HK CDF records.
%
% PROPOSAL: Separate function(s) for detecting sweeps, adding margin,
%           converting to science time(?).
%   * Detect sweeps using HK BDM for timestamps before specified limit.
%   * Detect sweeps using HK bias. Return value in HK time.
%   * Detect sweeps using L1R QUALITY_BITMASK. Return value in science time.
%   * Add time margins to time interval (t_before, t_after).
%   PROPOSAL: Replace file with class with static functions.
%     sweepdet
%
% PROPOSAL: Rename autodetect_sweeps() --> autodetect_sweeps_from_BIAS_HK



% The only BDM which sweeps use, but there may be other data too.
% I.e.
% BDM<>BDM_SWEEP_POSSIBLE ==> Not a sweep.
% BDM==BDM_SWEEP_POSSIBLE <== Sweep
BDM_SWEEP_POSSIBLE = 4;



% Time before which a sweep is equivalent to BDM=4. Inclusive threshold.
sbdaEndTt2000          = spdfcomputett2000(Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC'));
windowLengthPts        =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS');
% Minimum min-max difference for counting as sweep.
currentMmDiffMinimumTm =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MINMAX_DIFF_MINIMUM_TM');
windowMarginSec        =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MARGIN_SEC');



nCdfRecs = irf.assert.sizes(...
  hkTt2000,         [-1, 1], ...
  hkBdmFpa,         [-1, 1], ...
  hkBiasCurrentFpa, [-1, 3]);
% NOTE: Can not use integer MATLAB class (?) in settings.
assert(round(windowLengthPts) == windowLengthPts)
% NOTE: Must be at least two since SCDA requires comparison between
%       timestamps.
assert(windowLengthPts        >= 2)
assert(currentMmDiffMinimumTm >= 0)



hkBiasCurrent = hkBiasCurrentFpa.int2doubleNan();
bdm           = hkBdmFpa.int2doubleNan();
% Whether SBDA applies to (should be used for) records.
bSbdaApplies  = hkTt2000 <= sbdaEndTt2000;



%==========================
% Detect sweeps using SBDA
%==========================
isSweepingSbda = bSbdaApplies & (bdm == BDM_SWEEP_POSSIBLE);

%====================================
% Detect sweeps using sliding window
%====================================
isSweepingScda = false(size(isSweepingSbda));    % Preallocate
% NOTE: Will iterate zero times if window is longer than number of
%       records.
for i1 = 1:(nCdfRecs-(windowLengthPts-1))
  i2        = i1 + (windowLengthPts-1);
  iWindowAr = i1:i2;

  % NOTE: Treating all data in window combined, both in time AND over
  % channels/antennas.
  hkBiasCurrentWindow = hkBiasCurrent(i1:i2, :);

  hkBiasCurrentWindow(bdm(iWindowAr) ~= BDM_SWEEP_POSSIBLE, :) = NaN;
  minWindow = min(hkBiasCurrentWindow, [], 1);
  maxWindow = max(hkBiasCurrentWindow, [], 1);
  mmDiff    = maxWindow - minWindow;

  % NOTE: Can not reduce Only labelling
  if any(mmDiff >= currentMmDiffMinimumTm)
    isSweepingScda(iWindowAr) = isSweepingScda(iWindowAr) | ((bdm(iWindowAr) == BDM_SWEEP_POSSIBLE) & ~bSbdaApplies(iWindowAr));
  end
end

isSweeping           = isSweepingSbda | isSweepingScda;
isSweepingWithMargin = irf.utils.true_with_margin( ...
  hkTt2000, isSweeping, windowMarginSec * 1e9, windowMarginSec * 1e9);

isSweepingFpa        = bicas.utils.FPArray(isSweepingWithMargin);
end
