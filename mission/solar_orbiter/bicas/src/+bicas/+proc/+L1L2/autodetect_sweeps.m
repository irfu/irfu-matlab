% Try to autodetect sweeps.
%
% NOTE: This function should be temporary functionality while waiting
% for the long-term solution which is to use the relevant bit in L1/L1R
% QUALITY_BITMASK for detecting sweeps. Since it is a (hopefully)
% short-term solution, the functionality is also not as sophisticated
% (and presumably accurate) or configurable as it could be.
%
% NOTE: There is some unimportant "imprecision" in the window algorithm
% which can be seen as an unimportant bug. Records which are BDM<>4 (i.e.
% definately not sweep) within a window with currents exceeding the
% thresholds (for BDM=4) are labelled as sweeps.
%
%
% ALGORITHM
% =========
% PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC specifies a
% timestamps.
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
function isSweepingFpa = autodetect_sweeps(tt2000, bdmFpa, hkBiasCurrentFpa, Bso)
% TODO-DEC: Does having argument and return value FPAs make sense?
%           Should caller convert?
% PROPOSAL: Rename (and negate) BDM_SWEEP_POSSIBLE --> BDM_SWEEP_IMPOSSIBLE
%   CON: Unintuitive for time period when sweep can be deduced frmo BDM.
% PROPOSAL: Abbreviations for the two respective algorithms.
%   sweep
%   BDM, BDM=4
%   trick
%   autodetect
%   HK current, bias
%   algorithm, method
%   --
%   SBDA = Sweep(?) BDM Detection Algorithm  -- IMPLEMENTED
%   SCDA = Sweep Current Detection Algorithm -- IMPLEMENTED
%   SADA = Sweep AutoDetection Algorithm
%   BSDA = BDM Sweep Detection Algorithm
%     CON: Analogue "Current Sweep Detection Algorithm" is bad.

% The only BDM which sweeps use, but there may be other data too.
% I.e.
% BDM<>BDM_SWEEP_POSSIBLE ==> Not a sweep.
% BDM==BDM_SWEEP_POSSIBLE <== Sweep
BDM_SWEEP_POSSIBLE = 4;

% Time before which sweep <=> BDM=4.
sbdaEndTt2000          = spdfcomputett2000(Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC'));    % Inclusive threshold.
windowLengthPts        =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS');
currentMmDiffMinimumTm =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MINMAX_DIFF_MINIMUM_TM');   % Minimum value for counting as sweep.
windowMarginSec        =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MARGIN_SEC');

nCdfRecs = irf.assert.sizes(...
  tt2000,           [-1, 1], ...
  bdmFpa,           [-1, 1], ...
  hkBiasCurrentFpa, [-1, 3]);
% NOTE: Can not use integer in settings.
assert(round(windowLengthPts) == windowLengthPts)
% NOTE: Must be at least two since SCDA requires comparison between
%       timestamps.
assert(windowLengthPts        >= 2)
assert(currentMmDiffMinimumTm >= 0)

hkBiasCurrent = hkBiasCurrentFpa.int2doubleNan();
bdm           = bdmFpa.int2doubleNan();
% Whether SBDA applies to (should be used for) records.
bSbdaApplies  = tt2000 <= sbdaEndTt2000;

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
  i2 = i1 + (windowLengthPts-1);
  iWindowAr = i1:i2;

  % NOTE: Treating all data in window combined, both in time AND
  % over channels/antennas.
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
isSweepingWithMargin = irf.utils.true_with_margin(tt2000, isSweeping, windowMarginSec * 1e9);

isSweepingFpa        = bicas.utils.FPArray(isSweepingWithMargin);
end
